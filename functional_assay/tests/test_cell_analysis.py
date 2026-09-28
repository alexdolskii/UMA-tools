"""Scientific contracts, input validation, provenance, and run isolation."""

import contextlib
import csv
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
from functional_assay import calibration, cell_analysis, cell_imagej

from uma_tools import cli
from uma_tools.files import sha256_file


class OptionsTests(unittest.TestCase):
    def parse(self, *options):
        return cell_analysis.parse_args(["-i", "input.json", *options])

    def test_auto_and_all_manual_threshold_forms(self):
        for options, expected in (
            ([], None),
            (["-t"], (50, 65535)),
            (["-t", "100"], (100, 65535)),
            (["--threshold", "100", "5000"], (100, 5000)),
        ):
            with self.subTest(options=options):
                self.assertEqual(self.parse(*options).threshold, expected)
        self.assertEqual(self.parse().min_size_px, 5)

    def test_invalid_values_fail_before_any_analysis(self):
        for options in (
            ["-t", "nan"],
            ["-t", "inf"],
            ["-t", "-1"],
            ["-t", "65536"],
            ["-t", "60", "50"],
            ["-t", "1", "2", "3"],
            ["--min-size-px", "-1"],
            ["--min-size-um2", "nan"],
            ["--min-size-px", "inf"],
            ["--min-size-px", "5", "--min-size-um2", "20"],
        ):
            with self.subTest(options=options):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as error:
                        self.parse(*options)
                self.assertEqual(error.exception.code, 2)

    def test_pixel_and_physical_size_options_are_equivalent(self):
        pixel_area = 6.21480569402239**2
        manual = self.parse("--min-size-um2", str(5 * pixel_area))
        self.assertEqual(cell_analysis.minimum_sizes(manual, pixel_area)[0], 5)
        self.assertEqual(
            cell_analysis.minimum_sizes(self.parse(), 0.5 * 0.25),
            (5, 0.625),
        )

    def test_help_and_version_do_not_start_fiji(self):
        for option in ("--help", "--version"):
            script = (
                "import sys; from functional_assay.cell_analysis import main\n"
                f"try: main([{option!r}])\n"
                "except SystemExit as error: assert error.code == 0\n"
                "assert 'imagej' not in sys.modules\n"
                "assert 'scyjava' not in sys.modules\n"
            )
            result = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
                text=True,
                timeout=20,
            )
            self.assertEqual(result.returncode, 0, result.stderr)


class ParticleTests(unittest.TestCase):
    def test_area_removes_noise_and_keeps_holes_and_edge_particles(self):
        mask = np.zeros((20, 20), dtype=bool)
        mask[5:10, 5:10] = True
        mask[6:9, 6:9] = False  # ring: 16 pixels, not filled area of 25
        mask[0:2, 0:3] = True  # edge particle: 6 pixels
        mask[15, 15] = True  # noise: 1 pixel
        labels, sizes = cell_imagej.filter_particles(mask, 5)
        self.assertEqual(sorted(sizes.tolist()), [6, 16])
        self.assertEqual(np.count_nonzero(labels), 22)
        self.assertFalse(labels[7, 7])
        self.assertFalse(labels[15, 15])
        counted, areas = cell_imagej.filter_particles(
            labels != 0, 5, exclude_edges=True
        )
        self.assertEqual(areas.tolist(), [16])
        self.assertEqual(np.count_nonzero(counted), 16)

    def test_eight_connectivity_and_exact_size_boundary(self):
        mask = np.zeros((8, 8), dtype=bool)
        mask[2, 2] = mask[3, 3] = True
        labels, sizes = cell_imagej.filter_particles(mask, 2)
        self.assertEqual(sizes.tolist(), [2])
        self.assertEqual(labels[2, 2], labels[3, 3])
        _, sizes = cell_imagej.filter_particles(mask, 2.01)
        self.assertEqual(len(sizes), 0)

    def test_empty_mask_is_a_valid_zero_result(self):
        labels, sizes = cell_imagej.filter_particles(
            np.zeros((8, 8), dtype=bool), 5, exclude_edges=True
        )
        self.assertFalse(labels.any())
        self.assertEqual(len(sizes), 0)


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.files = [
            (index, Path(f"frame_{index}.nd2")) for index in range(9)
        ]

    def metadata(self, path):
        return {
            "filename": path.name,
            "width_px": 512,
            "height_px": 512,
            "slices": 3,
            "channels": 1,
            "frames": 1,
            "pixel_size_x_um": 0.5,
            "pixel_size_y_um": 0.25,
        }

    def test_all_nine_frames_are_checked_with_anisotropic_scale(self):
        with patch.object(
            calibration, "read_nd2_metadata", side_effect=self.metadata
        ) as read:
            result = calibration.read_well_calibration(self.files)
        self.assertEqual(read.call_count, 9)
        self.assertEqual(result["pixel_area_um2"], 0.125)
        self.assertEqual(
            [item["frame_index"] for item in result["tiles"]], list(range(9))
        )

    def test_duplicates_missing_and_inconsistent_metadata_are_rejected(self):
        for files in (self.files[:8], self.files[:8] + [self.files[0]]):
            with self.assertRaisesRegex(ValueError, "nine unique"):
                calibration.read_well_calibration(files)
        for key, bad in (
            ("pixel_size_x_um", 0.7),
            ("pixel_size_y_um", 0),
            ("pixel_size_x_um", float("nan")),
            ("slices", 4),
            ("channels", 2),
            ("frames", 2),
        ):

            def read(path):
                data = self.metadata(path)
                if path == self.files[-1][1]:
                    data[key] = bad
                return data

            with self.subTest(key=key, bad=bad):
                with patch.object(
                    calibration, "read_nd2_metadata", side_effect=read
                ):
                    with self.assertRaises(ValueError):
                        calibration.read_well_calibration(self.files)


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="UMA cells, ")
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.stitched = self.folder / "Stitched_Results"
        self.stitched.mkdir()
        self.path = self.stitched / "WellA1_stitched.tif"
        self.path.write_bytes(b"test input TIFF")
        self.files = []
        for index in range(9):
            raw = self.folder / f"x__WellA1_PointA1_{index:04d}_ChannelGFP.nd2"
            raw.write_bytes(b"original")
            self.files.append((index, raw))
        self.scale = {
            "pixel_size_x_um": 0.5,
            "pixel_size_y_um": 0.25,
            "pixel_area_um2": 0.125,
            "tiles": [{"filename": path.name} for _, path in self.files],
        }
        self.args = cell_analysis.parse_args(["-i", "input.json"])

    def metadata(self):
        return {
            "schema_version": 1,
            "overlap_percent": 30,
            "wells": {
                "WellA1": {
                    "status": "completed",
                    "sha256": sha256_file(self.path),
                    "width_px": 10,
                    "height_px": 8,
                    "calibration": self.scale,
                    "source_frames": self.scale["tiles"],
                }
            },
        }

    def result(self, empty=False):
        area = np.zeros((8, 10), dtype=bool)
        counted = area.copy()
        if not empty:
            area[2:4, 2:7] = True
            counted[2:4, 2:5] = True
        return {
            "area_mask": area,
            "counting_mask": counted,
            "object_areas_px2": np.array([] if empty else [6], dtype=int),
            "threshold_lower": 50,
            "threshold_upper": 65535,
            "width_px": 10,
            "height_px": 8,
        }

    def test_discovery_ignores_metadata_nonstitched_images_and_directories(
        self,
    ):
        for name in ("._WellA1_stitched.tif", "notes.txt", "projection.tif"):
            (self.stitched / name).write_bytes(b"ignored")
        (self.stitched / "WellA2_stitched.tif").mkdir()
        self.assertEqual(
            cell_analysis.discover_stitched(self.stitched),
            {"WellA1": self.path},
        )

    def test_provenance_binds_overlap_to_tiff_and_checks_scale(self):
        result = cell_analysis.verify_stitching_record(
            self.path, "WellA1", self.metadata(), self.scale
        )
        self.assertEqual(result["overlap_percent"], 30)
        legacy = cell_analysis.verify_stitching_record(
            self.path, "WellA1", None, self.scale
        )
        self.assertIsNone(legacy["overlap_percent"])
        data = self.metadata()
        self.path.write_bytes(b"different TIFF")
        with self.assertRaisesRegex(ValueError, "does not match"):
            cell_analysis.verify_stitching_record(
                self.path, "WellA1", data, self.scale
            )
        data = self.metadata()
        wrong = dict(self.scale, pixel_size_x_um=1)
        with self.assertRaisesRegex(ValueError, "calibration has changed"):
            cell_analysis.verify_stitching_record(
                self.path, "WellA1", data, wrong
            )

    def test_mask_area_and_counted_area_are_saved_separately_in_both_units(
        self,
    ):
        for empty in (False, True):
            with self.subTest(empty=empty):
                with (
                    patch.object(
                        cell_analysis,
                        "read_well_calibration",
                        return_value=self.scale,
                    ),
                    patch.object(
                        cell_imagej,
                        "analyze_image",
                        return_value=self.result(empty),
                    ),
                    patch.object(cell_imagej, "save_mask"),
                    patch.object(cell_imagej, "save_contours"),
                ):
                    row, _ = cell_analysis.measure_well(
                        "WellA1",
                        self.path,
                        self.files,
                        self.folder,
                        self.args,
                        None,
                    )
                self.assertEqual(row["Object_Count"], 0 if empty else 1)
                self.assertEqual(row["Mask_Area_px2"], 0 if empty else 10)
                self.assertEqual(row["Mask_Area_um2"], 0 if empty else 1.25)
                self.assertEqual(
                    row["Counted_Object_Area_um2"], 0 if empty else 0.75
                )
                self.assertEqual((row["Width_um"], row["Height_um"]), (5, 2))
                self.assertIsNone(row["Overlap_Percent"])
                with (self.folder / "WellA1_objects.csv").open() as stream:
                    objects = list(csv.DictReader(stream))
                self.assertEqual(len(objects), 0 if empty else 1)

    def test_run_keeps_old_results_and_reports_failed_well_without_fake_zero(
        self,
    ):
        (self.stitched / "WellB2_stitched.tif").write_bytes(b"incomplete well")
        metadata_path = self.stitched / "stitching_metadata.json"
        metadata_path.write_text(json.dumps(self.metadata()))
        outputs = []
        for _ in range(2):
            with (
                patch.object(
                    cell_analysis,
                    "read_well_calibration",
                    return_value=self.scale,
                ),
                patch.object(
                    cell_imagej, "analyze_image", return_value=self.result()
                ),
                patch.object(cell_imagej, "save_mask"),
                patch.object(cell_imagej, "save_contours"),
                contextlib.redirect_stdout(io.StringIO()),
            ):
                count, failures = cell_analysis.process_folder(
                    self.folder, self.args, Mock()
                )
            self.assertEqual((count, failures), (1, 1))
            outputs = sorted(self.folder.glob("Cell_Analysis_*"))
        self.assertEqual(len(outputs), 2)
        self.assertTrue(all(path.exists() for _, path in self.files))
        self.assertEqual(self.path.read_bytes(), b"test input TIFF")
        for output in outputs:
            status = json.loads((output / "run_status.json").read_text())
            self.assertEqual(status["status"], "partial")
            self.assertEqual(
                (output / metadata_path.name).read_bytes(),
                metadata_path.read_bytes(),
            )
            with (output / "Cell_Analysis_Summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["Overlap_Percent"], "30")
            self.assertEqual(rows[1]["Status"], "failed")
            self.assertEqual(rows[1]["Object_Count"], "")
            self.assertIn("nine unique", rows[1]["Error"])

    def test_imagej_failure_is_logged_in_new_output(self):
        with contextlib.redirect_stdout(io.StringIO()):
            count, failures = cell_analysis.process_folder(
                self.folder,
                self.args,
                Mock(side_effect=RuntimeError("JVM unavailable")),
            )
        self.assertEqual((count, failures), (0, 1))
        output = next(self.folder.glob("Cell_Analysis_*"))
        self.assertIn("JVM unavailable", (output / "run.log").read_text())
        self.assertEqual(
            json.loads((output / "run_status.json").read_text())["status"],
            "failed",
        )

    def test_explicit_appledouble_json_rejected_before_startup(self):
        args = cell_analysis.parse_args(
            ["-i", str(self.folder / "._input.json")]
        )
        with patch.object(cell_analysis, "initialize_imagej") as start:
            with patch.object(Path, "read_text") as read:
                with self.assertRaisesRegex(ValueError, "metadata"):
                    cell_analysis.run_analysis(args)
        start.assert_not_called()
        read.assert_not_called()

    def test_interrupted_run_is_never_marked_successful(self):
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(KeyboardInterrupt):
                cell_analysis.process_folder(
                    self.folder,
                    self.args,
                    Mock(side_effect=KeyboardInterrupt),
                )
        output = next(self.folder.glob("Cell_Analysis_*"))
        status = json.loads((output / "run_status.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["error"], "KeyboardInterrupt")

    def test_context_and_workers_close_on_success_and_failure(self):
        config = self.folder / "input.json"
        config.write_text(json.dumps({"folder_paths": [str(self.folder)]}))
        for fails in (False, True):
            context = Mock()

            def process(folder, args, ensure):
                ensure()
                return (0, 1) if fails else (1, 0)

            with (
                patch.object(
                    cell_analysis, "initialize_imagej", return_value=context
                ),
                patch.object(
                    cell_analysis, "process_folder", side_effect=process
                ),
                patch.object(cli, "_shutdown_imagej_workers") as workers,
                patch("scyjava.jvm_started", return_value=True),
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                exit_code = cell_analysis.main(["-i", str(config)])
            self.assertEqual(exit_code, 1 if fails else 0)
            context.dispose.assert_called_once_with()
            workers.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()
