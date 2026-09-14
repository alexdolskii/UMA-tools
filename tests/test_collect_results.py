"""End-to-end collection checks using temporary assay outputs, without Fiji."""

import contextlib
import csv
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from uma_tools import collect_results as collector


class CollectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "Исходные images"
        self.source.mkdir()
        self.config = self.root / "input paths.json"

    def config_for(self, *sources):
        self.config.write_text(
            json.dumps({"folder_paths": [str(path) for path in sources]}),
            encoding="utf-8",
        )
        return self.config

    def run_command(self, *sources):
        self.config_for(*(sources or (self.source,)))
        with contextlib.redirect_stdout(io.StringIO()) as terminal:
            code = collector.main(["-i", str(self.config)])
        return code, terminal.getvalue()

    def output(self, source=None):
        return sorted((source or self.source).glob("Combined_Results_*"))[-1]

    def status(self, source=None):
        return json.loads(
            (self.output(source) / "run_status.json").read_text(
                encoding="utf-8"
            )
        )

    def make_run(
        self,
        analysis,
        names=("field.one.nd2",),
        stamp="20260914_120000",
        source=None,
        angle="15",
        status="SUCCESS",
    ):
        source = source or self.source
        if analysis == "Alignment":
            run = source / f"Alignment_assay_results_angle_{angle}_{stamp}"
            path = run / "Analysis" / "Alignment_Summary.csv"
            columns = [
                "File_Name",
                "Number_of_Z_Stacks",
                "Z_Stack_Type",
                f"Percentage_Fibers_Aligned_Within_{angle}_Degree",
                "Orientation_Mode",
            ]
            # N/A is valid older alignment metadata, not a metric.
            rows = [
                [
                    Path(name).stem
                    + "_processed_orientation_distribution.csv",
                    "N/A",
                    "N/A",
                    65.5,
                    "Aligned",
                ]
                for name in names
            ]
        elif analysis == "Thickness":
            run = source / f"Thickness_assay_results_{stamp}"
            path = run / "Thickness_Summary.csv"
            columns = ["File_Name", "Area", "StdDev", "Min", "Max", "Median"]
            rows = [
                [
                    name,
                    1560.8484500537231,
                    0.900738580354334,
                    2.2360680103302,
                    7.0,
                    5.830951690673828,
                ]
                for name in names
            ]
        else:
            run = source / f"Area_assay_results_{stamp}_000123_42"
            path = run / "Fibronectin_Area_Summary.csv"
            columns = [
                "File_Name",
                "Image_ID",
                "FN_Positive_Pixels",
                "FN_Area_Percent",
                "FN_Area",
                "Threshold_Lower",
                "Threshold_Upper",
                "Projection_Method",
            ]
            rows = [
                [name, name, 12000, 30.5, 1200.55, 2000, "", "SUM"]
                for name in names
            ]
        path.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, lineterminator="\r\n")
        writer.writerow(columns)
        writer.writerows(rows)
        # Preserve BOM and CRLF exactly, without reserializing the CSV.
        path.write_bytes(b"\xef\xbb\xbf" + buffer.getvalue().encode("utf-8"))
        if analysis == "Area" and status is not None:
            (run / "run_status.json").write_text(
                json.dumps({"status": status}), encoding="utf-8"
            )
        return path

    def all_runs(
        self, names=("field.one.nd2",), source=None, stamp="20260914_120000"
    ):
        return {
            name: self.make_run(name, names, source=source, stamp=stamp)
            for name in ("Alignment", "Thickness", "Area")
        }

    def assert_no_summaries(self, output=None):
        output = output or self.output()
        self.assertEqual(list(output.glob("*_Summary.csv")), [])
        self.assertEqual(list(output.glob("*.partial")), [])

    def test_three_csvs_keep_bytes_and_exact_names_without_seq(self):
        names = ("клетка, first.field.nd2", "field_processed.tiff")
        paths = self.all_runs(names)
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        status = self.status()
        self.assertEqual(status["status"], "SUCCESS")
        self.assertEqual(status["image_check"], "MATCH")
        self.assertEqual(status["source_name"], self.source.name)
        self.assertIn(
            f"Combined_Results_{self.source.name}_", self.output().name
        )
        for original in paths.values():
            target = self.output() / f"{self.source.name}_{original.name}"
            self.assertEqual(target.read_bytes(), original.read_bytes())
        self.assertEqual(len(list(self.output().glob("*_Summary.csv"))), 3)
        with (self.output() / "image_check.csv").open(
            encoding="utf-8", newline=""
        ) as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual({row["Image_File_Name"] for row in rows}, set(names))
        self.assertTrue(all(row["Check"] == "MATCH" for row in rows))

    def test_newer_invalid_results_fall_back_independently(self):
        old = self.all_runs(stamp="20260914_100000")
        alignment = self.make_run("Alignment", stamp="20260914_140000")
        alignment.write_bytes(alignment.read_bytes().replace(b"65.5", b"nan"))
        self.make_run(
            "Thickness",
            names=("field.one.nd2", "field.one.nd2"),
            stamp="20260914_150000",
        )
        self.make_run("Area", stamp="20260914_160000", status="RUNNING")
        # Modifying/copying a historical CSV must not change the run ordering.
        for path in old.values():
            os.utime(path, (2000000000, 2000000000))
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        self.assertIn("nonfinite", terminal)
        self.assertIn("Duplicate image", terminal)
        self.assertIn("RUNNING", terminal)
        for analysis, path in old.items():
            self.assertEqual(
                self.status()["selected"][analysis]["summary"], str(path)
            )
        with (self.output() / "selection_report.csv").open(
            encoding="utf-8", newline=""
        ) as stream:
            records = list(csv.DictReader(stream))
        self.assertEqual(sum(row["Status"] == "INVALID" for row in records), 3)
        self.assertEqual(
            sum(row["Status"] == "SELECTED" for row in records), 3
        )

    def test_missing_two_one_and_zero_available_analyses(self):
        self.make_run("Alignment")
        self.make_run("Thickness")
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        self.assertIn("Found 2 of 3 analyses", terminal)
        self.assertIn(
            "Found 2 of 3 analyses",
            (self.output() / "run.log").read_text(encoding="utf-8"),
        )
        self.assertEqual(self.status()["missing_analyses"], ["Area"])
        self.assertEqual(len(list(self.output().glob("*_Summary.csv"))), 2)
        for count in (1, 0):
            with self.subTest(count=count):
                source = self.root / f"available_{count}"
                source.mkdir()
                if count:
                    self.make_run("Thickness", source=source)
                code, terminal = self.run_command(source)
                self.assertEqual(code, 0 if count else 1, terminal)
                self.assertEqual(
                    self.status(source)["image_check"],
                    "NOT_COMPARABLE" if count else "FAILED",
                )
                self.assertEqual(
                    len(list(self.output(source).glob("*_Summary.csv"))), count
                )

    def test_latest_mismatch_does_not_search_for_older_matching_combination(
        self,
    ):
        self.all_runs(("first.nd2",), stamp="20260914_100000")
        latest = self.make_run(
            "Thickness", ("second.nd2",), stamp="20260914_120000"
        )
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertEqual(
            self.status()["selected"]["Thickness"]["summary"], str(latest)
        )
        self.assertEqual(self.status()["status"], "VALIDATION_FAILED")
        self.assertIn("first.nd2", terminal)
        self.assertIn("second.nd2", terminal)
        self.assert_no_summaries()

    def test_equal_row_counts_with_different_extensions_fail(self):
        self.make_run("Thickness", ("field.nd2",))
        self.make_run("Area", ("field.tif",))
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("field.nd2", terminal)
        self.assertIn("field.tif", terminal)
        self.assert_no_summaries()

    def test_alignment_extension_ambiguity_is_not_silently_matched(self):
        self.all_runs(("field.nd2",))
        (self.source / "field.nd2").touch()
        (self.source / "field.tif").touch()
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("ambiguous", terminal)
        self.assertIn("field.tif", terminal)
        self.assert_no_summaries()

    def test_alignment_alone_checks_originals_and_unresolved_names(
        self,
    ):
        self.make_run("Alignment", ("sample_processed.v2.nd2",))
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("Cannot recover", terminal)
        self.assert_no_summaries()
        (self.source / "sample_processed.v2.nd2").touch()
        (self.source / "._sample_processed.v2.tif").touch()
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        self.assertEqual(self.status()["image_check"], "NOT_COMPARABLE")

    def test_empty_malformed_and_nonfinite_csvs_are_not_valid_results(self):
        replacements = [
            b"",
            b"File_Name,Area\nfield.nd2,1\n",
            b"File_Name,Area,StdDev,Min,Max,Median\n",
            b"File_Name,Area,StdDev,Min,Max,Median\nfield.nd2,1,2,3,4\n",
            b"File_Name,Area,StdDev,Min,Max,Median\nfield.nd2,1,2,3,4,5,6\n",
            b'File_Name,Area,StdDev,Min,Max,Median\n"unterminated',
            b"File_Name,Area,StdDev,Min,Max,Median\nfield.nd2,1,2,3,inf,5\n",
            b"File_Name,Area,StdDev,Min,Max,Median\nfield.nd2,1,2,3,4,\n",
            b"File_Name,Area,StdDev,Min,Max,Median,Median\nfield.nd2,1,2,3,4,5,6\n",
        ]
        for index, data in enumerate(replacements):
            with self.subTest(index=index):
                source = self.root / f"invalid_{index}"
                source.mkdir()
                path = self.make_run("Thickness", source=source)
                path.write_bytes(data)
                code, terminal = self.run_command(source)
                self.assertEqual(code, 1, terminal)
                self.assertEqual(self.status(source)["analyses_found"], 0)
                self.assert_no_summaries(self.output(source))

    def test_area_status_and_image_id_are_checked(self):
        for state in ("RUNNING", "ERROR", "CANCELLED", "VALIDATION_FAILED"):
            with self.subTest(state=state):
                source = self.root / state
                source.mkdir()
                self.make_run("Area", source=source, status=state)
                code, terminal = self.run_command(source)
                self.assertEqual(code, 1, terminal)
                self.assertIn(state, terminal)
        path = self.make_run("Area")
        path.write_bytes(
            path.read_bytes().replace(
                b"field.one.nd2,field.one.nd2", b"field.one.nd2,another.nd2"
            )
        )
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("Image_ID differs", terminal)

    def test_truncated_area_summary_disagrees_with_completion_counts(self):
        older = self.make_run("Area", stamp="20260914_100000")
        newer = self.make_run("Area", stamp="20260914_180000")
        (newer.parent / "run_status.json").write_text(
            json.dumps(
                {
                    "status": "SUCCESS",
                    "input_images": 2,
                    "processed_images": 2,
                    "generated_masks": 2,
                }
            ),
            encoding="utf-8",
        )
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        self.assertIn("row count", terminal)
        self.assertEqual(
            self.status()["selected"]["Area"]["summary"], str(older)
        )

    def test_same_latest_timestamp_is_reported_as_ambiguous(self):
        self.all_runs()
        self.make_run("Alignment", angle="10")
        code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("same latest timestamp", terminal)
        self.assert_no_summaries()

    def test_prior_collections_hidden_nested_and_partial_outputs_are_ignored(
        self,
    ):
        old = self.make_run("Thickness", stamp="20260914_100000")
        partial = self.make_run("Thickness", stamp="20260914_140000")
        partial.rename(partial.with_name("Thickness_Summary.partial.csv"))
        hidden = self.source / ".ignored"
        hidden.mkdir()
        self.all_runs(source=hidden, stamp="20260914_160000")
        nested = self.source / "unrelated"
        nested.mkdir()
        self.all_runs(source=nested, stamp="20260914_160000")
        linked = self.source / "Thickness_assay_results_20260914_180000"
        linked.symlink_to(old.parent, target_is_directory=True)
        hidden_run = self.source / "._Thickness_assay_results_20260914_190000"
        hidden_run.mkdir()
        (hidden_run / old.name).write_bytes(old.read_bytes())
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        first_output = self.output()
        before = {
            path.name: path.read_bytes()
            for path in first_output.iterdir()
            if path.is_file()
        }
        code, terminal = self.run_command()
        self.assertEqual(code, 0, terminal)
        self.assertEqual(len(list(self.source.glob("Combined_Results_*"))), 2)
        self.assertEqual(self.status()["analyses_found"], 1)
        self.assertEqual(
            self.status()["selected"]["Thickness"]["summary"], str(old)
        )
        self.assertEqual(
            before,
            {
                path.name: path.read_bytes()
                for path in first_output.iterdir()
                if path.is_file()
            },
        )

    def test_multiple_folders_continue_after_failure_and_deduplicate_paths(
        self,
    ):
        bad = self.root / "no analyses"
        bad.mkdir()
        second = self.root / "second originals"
        second.mkdir()
        self.all_runs()
        self.make_run("Thickness", source=second)
        code, terminal = self.run_command(
            bad, self.source, second, self.source / "."
        )
        self.assertEqual(code, 1, terminal)
        self.assertEqual(self.status(bad)["status"], "VALIDATION_FAILED")
        self.assertEqual(self.status()["status"], "SUCCESS")
        self.assertEqual(self.status(second)["source_name"], second.name)
        self.assertEqual(len(list(self.source.glob("Combined_Results_*"))), 1)
        self.assertIn("Skipping duplicate JSON folder", terminal)
        self.assertIn("2 folder(s) succeeded; 1 failed", terminal)

    def test_changed_source_is_not_published(self):
        paths = self.all_runs()
        original_publish = collector.publish_copies

        def change_then_publish(output, label, tables):
            paths["Thickness"].write_bytes(
                paths["Thickness"].read_bytes() + b"\n"
            )
            return original_publish(output, label, tables)

        with patch.object(
            collector, "publish_copies", side_effect=change_then_publish
        ):
            code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("Selected source changed", terminal)
        self.assert_no_summaries()

    def test_copy_failure_removes_already_published_summaries(self):
        self.all_runs()
        original_replace = Path.replace

        def fail_second_copy(path, target):
            if path.name.endswith("Thickness_Summary.csv.partial"):
                raise OSError("deliberate publication failure")
            return original_replace(path, target)

        with patch.object(Path, "replace", fail_second_copy):
            code, terminal = self.run_command()
        self.assertEqual(code, 1, terminal)
        self.assertIn("deliberate publication failure", terminal)
        self.assertEqual(self.status()["copied_csvs"], [])
        self.assertTrue((self.output() / "traceback.txt").is_file())
        self.assert_no_summaries()

    def test_metadata_json_rejected_before_read_and_startup_errors_logged(
        self,
    ):
        metadata = self.root / "._input.json"
        with patch.object(
            Path,
            "read_text",
            side_effect=AssertionError("Must not read metadata JSON"),
        ):
            with self.assertRaisesRegex(
                collector.ValidationError, "metadata JSON"
            ):
                collector.read_config(metadata)
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with contextlib.redirect_stdout(io.StringIO()):
                code = collector.main(["-i", str(metadata)])
            self.assertEqual(code, 1)
            outputs = list(
                self.root.glob("Combined_Results_configuration_error_*")
            )
            self.assertEqual(len(outputs), 1)
            self.assertIn(
                "metadata JSON", (outputs[0] / "run.log").read_text()
            )
            self.all_runs()
            missing = self.root / "nonexistent originals"
            code, terminal = self.run_command(missing, self.source)
            self.assertEqual(code, 1, terminal)
            self.assertFalse(missing.exists())
            self.assertEqual(self.status()["status"], "SUCCESS")
            self.assertTrue(
                list(
                    self.root.glob("Combined_Results_nonexistent originals_*")
                )
            )
        finally:
            os.chdir(previous)

    def test_installed_and_direct_commands_exit_without_imagej(self):
        self.all_runs()
        self.config_for(self.source)
        executable = Path(sys.executable).parent / "uma_collect_results"
        result = subprocess.run(
            [str(executable), "-i", str(self.config)],
            cwd=self.root,
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("3 analyses", result.stdout)
        self.assertIn("Collection finished", result.stdout)
        # -S disables site-packages; the launcher must use only stdlib.
        launcher = (
            Path(__file__).resolve().parents[1] / "code" / "collect_results.py"
        )
        result = subprocess.run(
            [sys.executable, "-S", str(launcher), "-i", str(self.config)],
            cwd=self.root,
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        script = (
            "import sys; from uma_tools import cli; "
            f"sys.argv=['uma_collect_results','-i',{str(self.config)!r}]; "
            "assert cli.collect_results() == 0; "
            "assert not {'imagej','jpype','scyjava','numpy','pandas'} "
            "& set(sys.modules)"
        )
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=self.root,
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.make_run("Thickness", ("different.nd2",), stamp="20260914_180000")
        result = subprocess.run(
            [str(executable), "-i", str(self.config)],
            cwd=self.root,
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("1 failed", result.stdout)


if __name__ == "__main__":
    unittest.main()
