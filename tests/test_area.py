"""Check source selection, error logs, and real SUM32 command results."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from uma_tools import area_analysis as area
from uma_tools.area_imagej import FLOAT32_MAX


class AreaInputTests(unittest.TestCase):
    def test_direct_inventory_needs_no_alignment_or_sequence_identifier(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in (
                "sample.nd2",
                "sample.tif",
                "another.TIFF",
                "._sample.nd2",
                ".hidden.tif",
                "notes.txt",
            ):
                (root / name).touch()
            nested = root / "previous_results.tiff"
            nested.mkdir()
            (nested / "old.tiff").touch()
            selected = area.original_inventory(root)
            self.assertEqual(
                [row["File_Name"] for row in selected],
                ["another.TIFF", "sample.nd2", "sample.tif"],
            )
            self.assertEqual(len({row["Image_ID"] for row in selected}), 3)
            run_id, output = area.new_output_folder(root)
            other_id, other_output = area.new_output_folder(root)
            self.assertNotEqual(output, other_output)
            self.assertNotEqual(run_id, other_id)
            self.assertEqual(len(area.original_inventory(root)), 3)

    def test_metadata_json_is_rejected_before_opening(self):
        args = area.parse_args(["-i", "._input_paths.json"])
        with patch.object(
            Path,
            "open",
            side_effect=AssertionError("Metadata JSON was opened"),
        ):
            with self.assertRaisesRegex(
                area.ValidationError, "macOS metadata"
            ):
                area.read_source_folders(args)

    def test_invalid_bounds_create_source_folder_diagnostics_without_java(
        self,
    ):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            manifest = root / "input.json"
            manifest.write_text(json.dumps({"folder_paths": [str(root)]}))
            with patch.object(area, "ImageJEngine") as engine:
                status = area.main(
                    [
                        "-i",
                        str(manifest),
                        "--channel",
                        "1",
                        "-t",
                        "5000",
                        "2000",
                    ]
                )
            self.assertEqual(status, 2)
            engine.assert_not_called()
            outputs = list(root.glob("Area_assay_results_*"))
            self.assertEqual(len(outputs), 1)
            self.assertTrue((outputs[0] / "errors.csv").is_file())
            self.assertTrue((outputs[0] / "traceback.txt").is_file())
            self.assertIn(
                "Upper threshold", (outputs[0] / "run.log").read_text()
            )
            self.assertEqual(
                json.loads((outputs[0] / "run_status.json").read_text())[
                    "status"
                ],
                "ERROR",
            )

    def test_bad_json_uses_results_folder_in_working_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            metadata = root / "._input.json"
            metadata.write_text("must not be read")
            command = Path(sys.executable).parent / "area_analysis"
            result = subprocess.run(
                [str(command), "-i", str(metadata)],
                cwd=root,
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 2, result.stderr)
            outputs = list(root.glob("Area_assay_results_*"))
            self.assertEqual(len(outputs), 1)
            self.assertIn(
                "macOS metadata", (outputs[0] / "run.log").read_text()
            )


@unittest.skipUnless(
    os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
    "Set UMA_RUN_IMAGEJ_TESTS=1 for real area command checks",
)
class AreaCommandImageJTests(unittest.TestCase):
    def write_multichannel_stack(self, path):
        import numpy as np
        import tifffile

        expected = np.array(
            [
                [0, 1999, 2000, 3000, 5000],
                [5001, 6000, 2000, 20, 2],
                [999, 1000, 10000, 4096, 60000],
                [65535, 70000, 80000, 100000, 131070],
            ],
            dtype=np.float32,
        )
        stack = np.zeros((2, 3, 4, 5), dtype=np.uint16)
        stack[:, 0] = 60000
        stack[0, 1] = expected // 2
        stack[1, 1] = expected - stack[0, 1]
        tifffile.imwrite(
            path,
            stack,
            imagej=True,
            resolution=(2, 4),
            metadata={"axes": "ZCYX", "spacing": 0.75, "unit": "um"},
        )
        return expected

    def run_command(self, root, manifest, bounds):
        command = Path(sys.executable).parent / "area_analysis"
        result = subprocess.run(
            [str(command), "-i", str(manifest), *bounds],
            cwd=root,
            input="2\n",
            capture_output=True,
            text=True,
            timeout=240,
        )
        self.assertEqual(
            result.stdout.count(
                "Enter fibronectin channel index (starting from 1):"
            ),
            1,
        )
        return result

    def test_sum_channel_bounds_calibration_and_reruns(self):
        import numpy as np
        import tifffile

        with tempfile.TemporaryDirectory(prefix="uma area ") as folder:
            root = Path(folder)
            source = root / "originals"
            source.mkdir()
            expected_sum = self.write_multichannel_stack(
                source / "sample.tiff"
            )
            (source / "._sample.tiff").write_bytes(b"Not an image")
            (source / "directory.tiff").mkdir()
            manifest = root / "input.json"
            manifest.write_text(json.dumps({"folder_paths": [str(source)]}))
            cases = [
                ((), 2000, FLOAT32_MAX),
                (("-t", "3000"), 3000, FLOAT32_MAX),
                (("--threshold", "2000", "5000"), 2000, 5000),
            ]
            for bounds, lower, upper in cases:
                with self.subTest(bounds=bounds):
                    before = set(source.glob("Area_assay_results_*"))
                    result = self.run_command(root, manifest, bounds)
                    self.assertEqual(
                        result.returncode, 0, result.stdout + result.stderr
                    )
                    outputs = set(source.glob("Area_assay_results_*")) - before
                    self.assertEqual(len(outputs), 1)
                    output = outputs.pop()
                    projection = next(
                        (output / "Projections_32bit").glob("*.tif")
                    )
                    values = tifffile.imread(projection)
                    self.assertEqual(values.dtype, np.dtype("float32"))
                    np.testing.assert_array_equal(values, expected_sum)
                    expected_mask = (expected_sum >= lower) & (
                        expected_sum <= upper
                    )
                    mask = tifffile.imread(
                        next((output / "Masks").glob("*.tif"))
                    )
                    np.testing.assert_array_equal(
                        mask, expected_mask.astype(np.uint8) * 255
                    )
                    with (output / "Fibronectin_Area_Summary.csv").open(
                        encoding="utf-8-sig"
                    ) as stream:
                        rows = list(csv.DictReader(stream))
                    self.assertEqual(len(rows), 1)
                    row = rows[0]
                    self.assertEqual(row["File_Name"], "sample.tiff")
                    self.assertEqual(int(row["Channel_Index"]), 2)
                    self.assertEqual(int(row["Number_of_Z_Stacks"]), 2)
                    self.assertEqual(
                        int(row["FN_Positive_Pixels"]),
                        int(expected_mask.sum()),
                    )
                    self.assertAlmostEqual(
                        float(row["FN_Area_Percent"]),
                        expected_mask.mean() * 100,
                    )
                    self.assertAlmostEqual(
                        float(row["FN_Area"]), expected_mask.sum() * 0.5 * 0.25
                    )
                    self.assertEqual(
                        float(row["Effective_Threshold_Upper"]), upper
                    )
                    status = json.loads(
                        (output / "run_status.json").read_text()
                    )
                    self.assertEqual(status["status"], "SUCCESS")
                    self.assertEqual(status["processed_images"], 1)
                    self.assertIn(
                        "ImageJ context and workers closed",
                        (output / "run.log").read_text(),
                    )
                    self.assertFalse(list(output.glob("*.partial.csv")))

    def test_failed_image_is_logged_and_next_folder_runs_before_exit_one(self):
        import numpy as np
        import tifffile

        with tempfile.TemporaryDirectory(prefix="uma area failure ") as folder:
            root = Path(folder)
            bad, good = root / "bad", root / "good"
            bad.mkdir()
            good.mkdir()
            tifffile.imwrite(
                bad / "single_channel.tif",
                np.ones((2, 4, 5), dtype=np.uint16),
                imagej=True,
                metadata={"axes": "ZYX"},
            )
            self.write_multichannel_stack(good / "valid.tiff")
            manifest = root / "input.json"
            manifest.write_text(
                json.dumps({"folder_paths": [str(bad), str(good)]})
            )
            result = self.run_command(root, manifest, ("-t", "2000", "5000"))
            self.assertEqual(
                result.returncode, 1, result.stdout + result.stderr
            )
            bad_output = next(bad.glob("Area_assay_results_*"))
            good_output = next(good.glob("Area_assay_results_*"))
            bad_status = json.loads(
                (bad_output / "run_status.json").read_text()
            )
            self.assertEqual(bad_status["status"], "VALIDATION_FAILED")
            self.assertEqual(bad_status["processed_images"], 0)
            self.assertEqual(bad_status["failed_images"], 1)
            self.assertFalse(
                (bad_output / "Fibronectin_Area_Summary.csv").exists()
            )
            self.assertTrue((bad_output / "errors.csv").is_file())
            self.assertTrue((bad_output / "traceback.txt").is_file())
            self.assertEqual(
                json.loads((good_output / "run_status.json").read_text())[
                    "status"
                ],
                "SUCCESS",
            )
            self.assertIn(
                "ImageJ context and workers closed",
                (bad_output / "run.log").read_text(),
            )


if __name__ == "__main__":
    unittest.main()
