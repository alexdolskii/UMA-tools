"""Cross-command recovery, new output layout and visible progress contracts."""

import argparse
import contextlib
import csv
import io
import json
import logging
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import openpyxl
from functional_assay import (
    cell_analysis,
    functional_report,
    report_data,
    stitching,
    survival_data,
    survival_report,
    workflow,
)
from test_functional_report import (
    create_analysis,
    create_template,
    example_rows,
)
from uma_tools.files import save_csv, save_json, sha256_file
from uma_tools.report_schema import ValidationError


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="UMA workflow ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source images"
        self.source.mkdir()
        self.config = self.root / "input.json"
        save_json(self.config, {"folder_paths": [str(self.source)]})
        self.stream = io.StringIO()
        redirect = contextlib.redirect_stdout(self.stream)
        redirect.__enter__()
        self.addCleanup(redirect.__exit__, None, None, None)

    def frames(self, well):
        for index in range(9):
            name = f"sample__{well}_PointA1_{index:04d}_ChannelGFP.nd2"
            (self.source / name).write_bytes(b"original frame")

    def stitched_partial(self):
        for well in ("WellA1", "WellA2", "WellA3"):
            self.frames(well)

        def fuse(folder, well, files, output, overlap, record):
            path = output / f"{well}_stitched.tif"
            path.write_bytes(well.encode())
            if well == "WellA2":
                raise RuntimeError("Synthetic fuse failure")
            record.update(
                output_file=path.name,
                width_px=64,
                height_px=64,
                slices=3,
                sha256=sha256_file(path),
            )

        with (
            patch.object(stitching, "stitch_well", side_effect=fuse) as mock,
            patch.object(
                stitching,
                "read_well_calibration",
                return_value={"pixel_size_x_um": 0.5, "pixel_size_y_um": 0.25},
            ),
        ):
            result = stitching.process_folder(self.source, 32.8, Mock())
        self.assertEqual(mock.call_count, 3)
        return result

    def test_stitch_failure_continues_and_partial_reaches_cell_count(self):
        stitched = self.stitched_partial()
        folder = Path(stitched["output"])
        self.assertEqual(stitched["status"], "PARTIAL")
        self.assertEqual(stitched["completed_wells"], 2)
        self.assertFalse((folder / "WellA2_stitched.tif").exists())
        self.assertTrue((folder / "WellA3_stitched.tif").is_file())
        self.assertIn("3/3 wells finished", self.stream.getvalue())

        def measure(well, path, *args):
            row = example_rows()[0]
            row.update(Well=well, File_Name=path.name)
            if well == "WellA1":
                for key in (
                    "Object_Count",
                    "Mask_Area_px2",
                    "Mask_Area_um2",
                    "Counted_Object_Area_px2",
                    "Counted_Object_Area_um2",
                ):
                    row[key] = 0
            return row, {}

        args = cell_analysis.parse_args(["-i", str(self.config)])
        with patch.object(
            cell_analysis, "measure_well", side_effect=measure
        ) as mock:
            self.assertEqual(
                cell_analysis.process_folder(self.source, args, Mock()), (2, 1)
            )
        self.assertEqual(
            [call.args[0] for call in mock.call_args_list],
            ["WellA1", "WellA3"],
        )
        selected = report_data.select_analysis(self.source, [])
        status = report_data.completed_status(selected / "run_status.json")
        rows = report_data.read_measurements(
            selected / "Cell_Analysis_Summary.csv", status
        )
        self.assertEqual(status["status"], "PARTIAL")
        self.assertEqual([row["Well"] for row in rows], ["A01", "A03"])
        self.assertEqual(rows[0]["Object_Count"], 0)
        excluded = read_csv(selected / "Processing_Exclusions.csv")
        self.assertEqual(excluded[0]["Well"], "WellA2")
        self.assertIn("Synthetic fuse failure", excluded[0]["Reason"])
        summary = selected / "Cell_Analysis_Summary.csv"
        summary.write_text(summary.read_text() + "\n")
        with self.assertRaisesRegex(ValidationError, "checksum"):
            report_data.read_measurements(summary, status)

    def test_empty_new_stitch_attempt_cannot_reuse_previous_tiffs(self):
        output = workflow.assay_directory(self.source) / "Stitched_Results"
        output.mkdir()
        (output / "WellA1_stitched.tif").write_bytes(b"stale")
        with patch.object(stitching, "stitch_well") as fuse:
            result = stitching.process_folder(self.source, 32.8)
        self.assertEqual(result["status"], "NO_INPUT")
        fuse.assert_not_called()
        self.assertFalse(list(output.glob("*_stitched.tif")))
        args = cell_analysis.parse_args(["-i", str(self.config)])
        start = Mock()
        self.assertEqual(
            cell_analysis.process_folder(self.source, args, start), (0, 1)
        )
        start.assert_not_called()
        analysis = next(output.parent.glob("Cell_Analysis_*"))
        self.assertEqual(
            json.loads((analysis / "run_status.json").read_text())["status"],
            "NO_INPUT",
        )

    def test_old_layout_is_never_selected_or_moved(self):
        analysis, _, _ = create_analysis(self.source)
        old = self.source / "Cell_Analysis_old_20270101_010101_000001"
        shutil.copytree(analysis, old)
        before = (old / "run_status.json").read_bytes()
        self.assertEqual(
            report_data.select_analysis(self.source, []), analysis
        )
        shutil.rmtree(analysis)
        with self.assertRaises(workflow.NoInputError):
            report_data.select_analysis(self.source, [])
        self.assertEqual((old / "run_status.json").read_bytes(), before)

    def test_failed_or_running_newer_runs_do_not_replace_finalized_partial(
        self,
    ):
        partial, _, _ = create_analysis(self.source, state="PARTIAL")
        for stamp, state in (
            ("20260929_010000_000001", "FAILED"),
            ("20260929_020000_000001", "RUNNING"),
        ):
            create_analysis(self.source, stamp, state)
        messages = []
        self.assertEqual(
            report_data.select_analysis(self.source, messages), partial
        )
        self.assertEqual(len(messages), 2)

    def test_partial_report_exports_successful_wells_and_audits_exclusions(
        self,
    ):
        analysis, _, _ = create_analysis(self.source, state="PARTIAL")
        create_template(analysis)
        args = argparse.Namespace(stats_unit="well", sheet=None)
        result = functional_report.process_folder(
            self.source, self.config, args
        )
        self.assertEqual(result["status"], "PARTIAL", result)
        output = Path(result["output"])
        self.assertEqual(output.parent, analysis.parent)
        self.assertEqual(len(read_csv(output / "Well_Data.csv")), 11)
        self.assertFalse(list(analysis.glob("Functional_Report_*")))
        excluded = read_csv(output / "Processing_Exclusions.csv")
        self.assertEqual(len(excluded), 1)
        self.assertEqual(excluded[0]["Well"], "WellC07")
        plots = json.loads((output / "plot_manifest.json").read_text())
        self.assertTrue(all(len(plot["Wells"]) == 11 for plot in plots))
        statistics = read_csv(output / "Statistics.csv")
        self.assertTrue(any(row["Treatment_N"] == "2" for row in statistics))
        book = openpyxl.load_workbook(output / result["workbook"])
        try:
            self.assertEqual(book["Processing Exclusions"].max_row, 2)
            self.assertIn(
                ("Report status", "PARTIAL"), list(book["Run Details"].values)
            )
        finally:
            book.close()

    def test_failed_wells_cannot_hide_measurements_or_unregistered_errors(
        self,
    ):
        analysis, rows, status = create_analysis(self.source, state="PARTIAL")
        path = analysis / "Cell_Analysis_Summary.csv"
        rows[-1]["Object_Count"] = 0
        save_csv(path, cell_analysis.SUMMARY_COLUMNS, rows)
        with self.assertRaisesRegex(ValidationError, "failed well contains"):
            report_data.read_measurements(path, status)
        rows[-1].pop("Object_Count")
        rows[-1]["Error"] = "unrecorded different error"
        save_csv(path, cell_analysis.SUMMARY_COLUMNS, rows)
        with self.assertRaises(ValidationError):
            report_data.read_measurements(path, status)

    def test_rotating_journals_keep_tracebacks_and_do_not_mix_commands(self):
        parent = workflow.assay_directory(self.source)
        for name in ("first", "second", "third"):
            (parent / name).mkdir()
        log = workflow.RunLog(parent / "first", self.source, "cell_count")
        try:
            try:
                raise RuntimeError("Diagnostic marker from run one")
            except RuntimeError as error:
                log.record_error("WellA01", error)
            logging.warning("External library warning")
        finally:
            log.close()
        journal = parent / "UMA_Logs" / "2_cell_count.log"
        before = journal.read_text()
        log = workflow.RunLog(parent / "second", self.source, "cell_count")
        log.event("SUCCESS", "Run two", "New marker")
        log.close()
        self.assertNotIn("run one", journal.read_text())
        archive = list((journal.parent / "archive").glob("2_cell_count_*.log"))
        self.assertEqual(len(archive), 1)
        self.assertEqual(archive[0].read_text(), before)
        self.assertIn("Traceback", before)
        self.assertIn("External library warning", before)
        saved = journal.read_bytes()
        log = workflow.RunLog(
            parent / "third", self.source, "functional_report"
        )
        log.close()
        self.assertEqual(journal.read_bytes(), saved)
        self.assertNotIn("\x1b", self.stream.getvalue())

    def test_live_progress_keeps_well_count_during_operation_and_stops_thread(
        self,
    ):
        class Terminal(io.StringIO):
            def isatty(self):
                return True

        output = workflow.assay_directory(self.source) / "run"
        output.mkdir()
        terminal = Terminal()
        with contextlib.redirect_stdout(terminal):
            log = workflow.RunLog(output, self.source, "stitching")
            log.phase(2, 3, "Stitch", finished=1, count=3)
            workflow.activity("WellA02: fuse nine stacks")
            self.assertIn("1/3 wells finished", log.progress.label)
            self.assertIn("WellA02", log.progress.label)
            self.assertTrue(log.progress.thread.is_alive())
            log.close()
            self.assertFalse(log.progress.thread.is_alive())
        self.assertIn("elapsed", terminal.getvalue())

    def test_folder_io_failure_continues_next_folder_and_deduplicates(self):
        second = self.root / "second"
        second.mkdir()
        save_json(
            self.config,
            {"folder_paths": [str(self.source), str(second), str(second)]},
        )
        args = cell_analysis.parse_args(["-i", str(self.config)])
        with (
            patch.object(
                cell_analysis,
                "process_folder",
                side_effect=[PermissionError("Read-only output"), (1, 0)],
            ) as process,
            contextlib.redirect_stderr(io.StringIO()),
        ):
            result = cell_analysis.run_analysis(args)
        self.assertEqual(result, 1)
        self.assertEqual(
            [call.args[0] for call in process.call_args_list],
            [self.source, second],
        )
        journal = (
            self.source
            / "uma_functional_assay"
            / "UMA_Logs"
            / "2_cell_count.log"
        )
        self.assertIn("Read-only output", journal.read_text())

    def test_ctrl_c_returns_130_and_records_cancelled_stitch(self):
        self.frames("WellA1")
        with (
            patch.object(
                stitching, "initialize_imagej", side_effect=KeyboardInterrupt
            ),
            patch("uma_tools.cli._shutdown_imagej_workers") as shutdown,
            contextlib.redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(stitching.main(["-i", str(self.config)]), 130)
        shutdown.assert_called_once()
        status = json.loads(
            (
                self.source
                / "uma_functional_assay"
                / "Stitched_Results"
                / "run_status.json"
            ).read_text()
        )
        self.assertEqual(status["status"], "CANCELLED")


class SurvivalRecoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="UMA day recovery ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.template = create_template(self.root)
        self.path = self.root / "survival.json"
        points = []
        self.analyses = {}
        for day in (1, 3, 5):
            source = self.root / f"images_{day}"
            source.mkdir()
            analysis, _, _ = create_analysis(source)
            self.analyses[day] = analysis
            points.append({"day": day, "folder": str(source)})
        save_json(
            self.path,
            {
                "experiment_name": "Day recovery",
                "plate_template": str(self.template),
                "output_dir": str(self.root),
                "baseline_day": 1,
                "difference_days": [3, 5],
                "timepoints": points,
            },
        )

    def report(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return survival_report.run_report(
                survival_data.read_config(self.path), self.path, True
            )

    def test_whole_later_day_missing_preserves_other_day_pairs_without_zeros(
        self,
    ):
        shutil.rmtree(self.analyses[3].parent.parent)
        result = self.report()
        self.assertEqual(result["status"], "PARTIAL", result)
        self.assertEqual(result["missing_days"], [3])
        output = Path(result["output"])
        self.assertEqual(output.parent, self.root / "uma_functional_assay")
        rows = read_csv(output / "Changes_by_Well.csv")
        missing = [row for row in rows if row["Day"] == "3"]
        self.assertEqual(len(missing), 24)
        self.assertTrue(
            all(
                row["Delta"] == "" and row["Status"] == "MISSING_DAY"
                for row in missing
            )
        )
        remaining = [row for row in rows if row["Day"] == "5"]
        self.assertTrue(
            all(
                row["Status"] == "PAIRED" and float(row["Delta"]) == 0
                for row in remaining
            )
        )
        statistics = read_csv(output / "Statistics.csv")
        self.assertEqual({row["Family_Size"] for row in statistics}, {"4"})
        self.assertTrue(
            all(
                row["Status"] == "Not tested"
                for row in statistics
                if row["Day"] == "3"
            )
        )
        self.assertEqual(result["paired_well_changes"], 12)
        self.assertEqual(len(result["plots"]), 6)

    def test_partial_day_uses_latest_successes_and_reduces_only_affected_pairs(
        self,
    ):
        source = self.analyses[3].parent.parent
        latest, _, _ = create_analysis(
            source, "20260930_010000_000001", "PARTIAL"
        )
        result = self.report()
        self.assertEqual(result["status"], "PARTIAL", result)
        self.assertEqual(result["paired_well_changes"], 23)
        selected = next(row for row in result["selections"] if row["Day"] == 3)
        self.assertEqual(selected["Selected_Analysis"], str(latest))
        self.assertEqual(selected["Completed_Wells"], 11)
        excluded = read_csv(
            Path(result["output"]) / "Processing_Exclusions.csv"
        )
        self.assertEqual(
            [(row["Day"], row["Well"]) for row in excluded], [("3", "WellC07")]
        )

    def test_absent_baseline_stops_changes_and_saves_day_diagnostics(self):
        shutil.rmtree(self.analyses[1].parent.parent)
        result = self.report()
        self.assertEqual(result["status"], "FAILED")
        self.assertIn("Baseline day 1", result["error"])
        output = Path(result["output"])
        self.assertFalse(list(output.glob("*.png")))
        selections = read_csv(output / "Selected_Analyses.csv")
        self.assertEqual(len(selections), 3)
        self.assertEqual(selections[0]["Status"], "NO_INPUT")

    def test_partial_baseline_keeps_only_wells_available_on_both_days(self):
        source = self.analyses[1].parent.parent
        create_analysis(source, "20260930_010000_000001", "PARTIAL")
        result = self.report()
        self.assertEqual(result["status"], "PARTIAL", result)
        self.assertEqual(result["paired_well_changes"], 22)
        self.assertEqual(result["missing_days"], [])
        rows = read_csv(Path(result["output"]) / "Changes_by_Well.csv")
        affected = [row for row in rows if row["Well"] == "C07"]
        self.assertEqual(len(affected), 4)
        self.assertTrue(
            all(
                row["Status"] == "MISSING_BASELINE" and row["Delta"] == ""
                for row in affected
            )
        )

    def test_baseline_without_any_mapped_well_stops_changes(self):
        analysis = self.analyses[1]
        status = json.loads((analysis / "run_status.json").read_text())
        rows = example_rows()[:1]
        rows[0].update(Well="WellH12", File_Name="WellH12_stitched.tif")
        status.update(
            completed_wells=1, wells={"WellH12": {"status": "completed"}}
        )
        save_csv(
            analysis / "Cell_Analysis_Summary.csv",
            cell_analysis.SUMMARY_COLUMNS,
            rows,
        )
        save_json(analysis / "run_status.json", status)
        result = self.report()
        self.assertEqual(result["status"], "FAILED")
        self.assertIn("no usable, mapped wells", result["error"])


if __name__ == "__main__":
    unittest.main()
