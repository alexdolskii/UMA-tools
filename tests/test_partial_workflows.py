"""Recoverable failures, auditable collection and real report handoff."""

import contextlib
import csv
import io
import json
import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import openpyxl
import test_collect_results as collection_fixtures
import test_report as report_fixtures
import tifffile

from uma_tools import alignment_analysis as alignment
from uma_tools import area_analysis as area
from uma_tools import cli, collect_results
from uma_tools import thickness_analysis as thickness
from uma_tools.image_run import ImageRun
from uma_tools.progress import (
    Cancelled,
    CommandSession,
    CompactProgress,
    ask_channel,
    confirm_start,
    folder_scope,
    phase,
)
from uma_tools.run import RunLog


def load_status(path):
    return json.loads((path / "run_status.json").read_text())


def audited_summary(path, names, failures):
    run_dir = (
        path.parent.parent if path.parent.name == "Analysis" else path.parent
    )
    audit = ImageRun(
        run_dir, run_dir.parent.parent, names, {"threshold": 2500}
    )
    for name in names:
        with audit.attempt(name, "Measurement", final=True):
            if name in failures:
                raise ValueError("Deliberate unreadable image")
    audit.finish(path)


class PartialCollectionTests(unittest.TestCase):
    def setUp(self):
        self.fixture = collection_fixtures.CollectionTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.names = (
            "sample_WellB02_1.nd2",
            "sample_WellB03_1.nd2",
            "broken_WellC02_1.nd2",
        )
        for name in self.names:
            (self.fixture.source / name).touch()
        self.paths = {
            assay.name: self.fixture.make_run(assay.name, self.names)
            for assay in collect_results.ASSAYS
        }

    def add_partial(self):
        path = self.fixture.make_run(
            "Area",
            self.names[:2],
            stamp="20261001_120000",
            threshold_label="2500",
        )
        audited_summary(path, self.names, {self.names[2]})
        return path

    def test_latest_partial_collects_matched_copies_and_reports_exclusions(
        self,
    ):
        partial = self.add_partial()
        original_bytes = {
            name: path.read_bytes() for name, path in self.paths.items()
        }
        with contextlib.redirect_stdout(io.StringIO()):
            code, _ = self.fixture.run_command()
        self.assertEqual(code, 0)
        combined = self.fixture.output()
        state = self.fixture.status()
        self.assertEqual(state["selected"]["Area"]["summary"], str(partial))
        self.assertEqual(state["excluded_images"], 1)
        self.assertEqual(state["retained_images"], 2)
        for item in state["copied_csvs"]:
            with Path(item["path"]).open(
                encoding="utf-8-sig", newline=""
            ) as stream:
                self.assertEqual(len(list(csv.DictReader(stream))), 2)
        for name, path in self.paths.items():
            self.assertEqual(path.read_bytes(), original_bytes[name])
        report_fixtures.ReportFixture.make_template(combined / "my plate.xlsx")
        with contextlib.redirect_stdout(io.StringIO()):
            result = cli.report.__wrapped__(
                ["-i", str(self.fixture.config), "--fn-threshold", "20"]
            )
        self.assertEqual(result, 0)
        output = next(combined.parent.glob("UMA_Report_*"))
        report_status = load_status(output)
        self.assertEqual(report_status["total_images"], 2)
        self.assertEqual(report_status["processing_excluded_images"], 1)
        self.assertEqual(report_status["generated_plots"], 6)
        workbook = openpyxl.load_workbook(
            report_status["workbook"], read_only=True
        )
        try:
            rows = list(workbook["Processing Exclusions"].values)
            self.assertEqual(rows[1][0], self.names[2])
            self.assertIn("Deliberate unreadable image", rows[1][1])
            self.assertEqual(workbook["Merged Data"].max_row, 3)
        finally:
            workbook.close()
        text = (
            self.fixture.source / "uma_assay/UMA_Logs/5_report.log"
        ).read_text()
        self.assertIn(self.names[2], text)
        self.assertIn("Processing exclusions", text)

    def test_unexplained_missing_image_fails_without_matching_fallback(
        self,
    ):
        self.add_partial()
        self.fixture.make_run(
            "Thickness", self.names[:1], stamp="20261001_130000"
        )
        code, terminal = self.fixture.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("Unexplained missing row", terminal)
        self.assertEqual(self.fixture.status()["copied_csvs"], [])

    def test_partial_requires_complete_audit_and_matching_summary_fingerprint(
        self,
    ):
        partial = self.add_partial()
        status_path = partial.parent / "run_status.json"
        valid = json.loads(status_path.read_text())
        for mutation in ("pending", "fingerprint", "unregistered"):
            with self.subTest(mutation=mutation):
                state = json.loads(json.dumps(valid))
                if mutation == "pending":
                    state["images"][-1]["Status"] = "PENDING"
                elif mutation == "fingerprint":
                    state["summary_sha256"] = "0" * 64
                else:
                    del state["image_status_schema"]
                status_path.write_text(json.dumps(state))
                code, _ = self.fixture.run_command()
                self.assertEqual(code, 0)
                selected = self.fixture.status()["selected"]["Area"]["summary"]
                self.assertEqual(selected, str(self.paths["Area"]))

    def test_report_rejects_tampered_exclusion_table(self):
        self.add_partial()
        self.fixture.run_command()
        combined = self.fixture.output()
        report_fixtures.ReportFixture.make_template(combined / "plate.xlsx")
        (combined / "processing_exclusions.csv").write_text("changed")
        with contextlib.redirect_stdout(io.StringIO()):
            result = cli.report.__wrapped__(["-i", str(self.fixture.config)])
        self.assertEqual(result, 1)
        status = load_status(next(combined.parent.glob("UMA_Report_*")))
        self.assertIn("exclusions", status["error"])
        self.assertNotIn("workbook", status)


class ImageFailureTests(unittest.TestCase):
    def test_thickness_keeps_good_images_after_exception_and_none_result(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            names = ("bad.tiff", "empty.tiff", "good.tiff")
            for name in names:
                (source / name).touch()
            result = dict(
                File_Name="good.tiff",
                Area=10,
                StdDev=1,
                Min=1,
                Max=4,
                Median=2,
            )

            def process(*args):
                if args[6] == "bad.tiff":
                    raise ValueError("Broken image")
                return result if args[6] == "good.tiff" else None

            with patch.object(
                thickness, "process_single_file", side_effect=process
            ) as calculate:
                state = thickness.process_single_folder(
                    None, None, None, None, None, str(source), ".tiff", 1
                )
            self.assertEqual(state, "PARTIAL")
            self.assertEqual(calculate.call_count, 3)
            status = load_status(
                next(source.glob("uma_assay/Thickness_assay_results_*"))
            )
            self.assertEqual(status["processed_images"], 1)
            self.assertEqual(status["failed_images"], 2)
            self.assertEqual(status["unprocessed_images"], 0)

    def test_alignment_failed_orientation_is_excluded_from_summary(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            for name in ("bad.nd2", "good.nd2"):
                (source / name).touch()

            def prepare(folder, output, *_):
                y, x = np.mgrid[:32, :32]
                pixels = (
                    200
                    * np.exp(-(((x - 15) / 3) ** 2))
                    * np.exp(-(((y - 15) / 12) ** 2))
                ).astype(np.uint8)
                tifffile.imwrite(Path(output) / "good_processed.tif", pixels)
                (Path(output) / "bad_processed.tif").write_bytes(b"broken")
                return {
                    f"{stem}_processed": {
                        "original_filename": f"{stem}.nd2",
                        "number_of_z_stacks": 3,
                        "z_stack_type": "slices",
                    }
                    for stem in ("bad", "good")
                }

            with patch.object(alignment, "process_part1", side_effect=prepare):
                state = alignment.process_folder(
                    str(source), 1, 15, 32, 32, None
                )
            self.assertEqual(state, "PARTIAL")
            output = next(source.glob("uma_assay/Alignment_assay_results_*"))
            status = load_status(output)
            self.assertEqual(status["processed_images"], 1)
            with (output / "Analysis/Alignment_Summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 1)
            self.assertTrue(rows[0]["File_Name"].startswith("good_"))

    def test_all_commands_mark_empty_sources_no_input_without_starting_java(
        self,
    ):
        for entry, module, answers, extra in (
            (cli.alignment.__wrapped__, alignment, ["1", "y"], []),
            (cli.thickness.__wrapped__, thickness, ["2", "1", "y"], []),
            (cli.area.__wrapped__, area, [], ["--channel", "1"]),
        ):
            with (
                self.subTest(command=entry.__name__),
                tempfile.TemporaryDirectory() as temporary,
            ):
                source = Path(temporary)
                config = source / "input.json"
                config.write_text(json.dumps({"folder_paths": [str(source)]}))
                target = (
                    "ImageJEngine" if module is area else "initialize_imagej"
                )
                with (
                    patch.object(module, target) as initialize,
                    patch("builtins.input", side_effect=answers),
                ):
                    code = entry(["-i", str(config), *extra])
                self.assertEqual(code, 1)
                initialize.assert_not_called()
                output = next(source.glob("uma_assay/*assay_results_*"))
                self.assertEqual(load_status(output)["status"], "NO_INPUT")

    def test_folder_failure_does_not_stop_next_folder_and_disposes_context(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            folders = [source / "one", source / "two"]
            for folder in folders:
                folder.mkdir()
                (folder / "image.tiff").touch()
            config = source / "input.json"
            config.write_text(
                json.dumps({"folder_paths": list(map(str, folders))})
            )
            gateway = Mock()
            with (
                patch.object(
                    alignment, "initialize_imagej", return_value=gateway
                ),
                patch.object(
                    alignment,
                    "process_folder",
                    side_effect=[ValueError("first failed"), "SUCCESS"],
                ) as process,
                patch("builtins.input", side_effect=["1", "y"]),
            ):
                code = cli.alignment.__wrapped__(["-i", str(config)])
            self.assertEqual(code, 1)
            self.assertEqual(process.call_count, 2)
            gateway.dispose.assert_called_once_with()
            first = (
                folders[0] / "uma_assay/UMA_Logs/1_alignment.log"
            ).read_text()
            second = (
                folders[1] / "uma_assay/UMA_Logs/1_alignment.log"
            ).read_text()
            self.assertIn("first failed", first)
            self.assertNotIn("first failed", second)


class JournalTests(unittest.TestCase):
    def test_rotation_routing_restoration_and_phase_capture(self):
        original_handlers = logging.getLogger().handlers[:]
        original_level = logging.getLogger().level
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            a, b = root / "a", root / "b"
            a.mkdir()
            b.mkdir()
            for index in range(2):
                with CommandSession("area", [a, b], "input.json") as session:
                    with folder_scope(a):
                        phase(f"a phase {index}")
                        out = a / f"run{index}"
                        out.mkdir()
                        log = RunLog(out)
                        try:
                            log.event("WARNING", "Measure", "a warning")
                        finally:
                            log.close()
                    with folder_scope(b):
                        phase(f"b phase {index}")
                    session.exit_code = 0
            latest = (a / "uma_assay/UMA_Logs/3_area.log").read_text()
            archived = list(
                (a / "uma_assay/UMA_Logs/archive").glob("3_area_*.log")
            )
            self.assertEqual(len(archived), 1)
            self.assertIn("a phase 0", archived[0].read_text())
            self.assertNotIn("a phase 0", latest)
            self.assertIn("a phase 1", latest)
            self.assertNotIn("b phase", latest)
            self.assertEqual(latest.count("a warning"), 1)
        self.assertEqual(logging.getLogger().handlers, original_handlers)
        self.assertEqual(logging.getLogger().level, original_level)

    def test_prompts_retry_and_cancel_without_java(self):
        with patch("builtins.input", side_effect=["wrong", "0", "-2", "3"]):
            self.assertEqual(ask_channel(), 3)
        with patch("builtins.input", side_effect=["wrong", "n"]):
            with self.assertRaises(Cancelled):
                confirm_start()
        for answer in ("q", EOFError(), KeyboardInterrupt()):
            with tempfile.TemporaryDirectory() as temporary:
                source = Path(temporary)
                config = source / "input.json"
                config.write_text(json.dumps({"folder_paths": [str(source)]}))
                with (
                    patch("builtins.input", side_effect=[answer]),
                    patch.object(alignment, "initialize_imagej") as initialize,
                ):
                    code = cli.alignment.__wrapped__(["-i", str(config)])
                self.assertEqual(code, 130)
                initialize.assert_not_called()
                log = (
                    source / "uma_assay/UMA_Logs/1_alignment.log"
                ).read_text()
                self.assertIn("CANCELLED", log)
                self.assertNotIn("Traceback", log)

    def test_redirected_progress_contains_no_ansi_or_fake_percentage(self):
        stream = io.StringIO()
        with CompactProgress(stream) as progress:
            progress.update("Reading image")
            progress.message("WARNING: a damaged image")
        self.assertEqual(stream.getvalue(), "WARNING: a damaged image\n")
        self.assertIsNone(progress.thread)

    def test_area_shutdown_is_recorded_once_in_each_source_journal(self):
        with tempfile.TemporaryDirectory() as temporary:
            sources = [Path(temporary) / name for name in ("one", "two")]
            outputs = []
            for source in sources:
                output = source / "uma_assay" / "Area_assay_results_test"
                output.mkdir(parents=True)
                (output / "run_status.json").write_text(
                    json.dumps(
                        {"source_folder": str(source), "status": "SUCCESS"}
                    )
                )
                outputs.append(output)
            with CommandSession("area", sources, "input.json") as session:
                with patch.object(area, "shutdown_imagej_workers"):
                    self.assertTrue(area.finish_imagej([None], outputs))
                session.exit_code = 0
            for source in sources:
                log = (source / "uma_assay/UMA_Logs/3_area.log").read_text()
                self.assertEqual(
                    log.count("ImageJ context and workers closed"), 1
                )

    def test_unwritable_journal_skips_only_that_folder(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bad, good = root / "bad", root / "good"
            bad.mkdir()
            good.mkdir()
            (bad / "uma_assay").mkdir()
            (bad / "uma_assay/UMA_Logs").write_text("not a directory")
            with CommandSession(
                "report", [bad, good], "input.json"
            ) as session:
                with self.assertRaises(OSError):
                    with folder_scope(bad):
                        self.fail("A folder without its journal was processed")
                with folder_scope(good):
                    phase("Continue good folder")
                session.exit_code = 1
            self.assertIn(
                "Continue good folder",
                (good / "uma_assay/UMA_Logs/5_report.log").read_text(),
            )
