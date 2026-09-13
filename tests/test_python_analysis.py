"""Exercise Python processing, input filtering, and cleanup without starting Java."""

from contextlib import redirect_stdout
import io
import json
import logging
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import tifffile

from uma_tools import alignment_analysis as alignment
from uma_tools import thickness_analysis as thickness


class PythonAnalysisTests(unittest.TestCase):
    def test_hidden_files_are_not_counted(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            images = root / "images"
            images.mkdir()
            (images / "sample.tif").touch()
            (images / "._sample.tif").touch()
            (images / ".hidden.tif").touch()
            (images / "subfolder").mkdir()
            manifest = root / "input.json"
            manifest.write_text(json.dumps({"folder_paths": [str(images)]}))
            for module in (alignment, thickness):
                output = io.StringIO()
                with redirect_stdout(output):
                    self.assertEqual(module.get_folder_paths(str(manifest)), [str(images)])
                self.assertIn("Number of files: 1", output.getvalue())
                with self.assertRaises(ValueError):
                    module.get_folder_paths(str(root / "._input.json"))

    def test_missing_input_does_not_start_java(self):
        for module, function in ((alignment, alignment.main_fibronectin_processing),
                                 (thickness, thickness.main)):
            with patch.object(module, "initialize_imagej") as initialize:
                with self.assertRaises(FileNotFoundError):
                    function("/nonexistent/uma-input.json")
                initialize.assert_not_called()

    def test_alignment_disposes_imagej_after_error(self):
        gateway = Mock()
        with patch.object(alignment, "get_folder_paths", return_value=["images"]), \
             patch.object(alignment, "initialize_imagej", return_value=gateway), \
             patch("builtins.input", side_effect=["1", "y"]), \
             patch.object(alignment, "process_folder", side_effect=RuntimeError("test failure")):
            with self.assertRaisesRegex(RuntimeError, "test failure"):
                alignment.main_fibronectin_processing("input.json")
        gateway.dispose.assert_called_once_with()

    def test_thickness_disposes_imagej_after_error(self):
        gateway = Mock()
        with patch.object(thickness, "get_folder_paths", return_value=["images"]), \
             patch.object(thickness, "initialize_imagej", return_value=gateway), \
             patch.object(thickness, "import_java_classes", return_value=(Mock(),) * 6), \
             patch("builtins.input", side_effect=["2", "1", "y"]), \
             patch.object(thickness, "process_all_folders", side_effect=RuntimeError("test failure")):
            with self.assertRaisesRegex(RuntimeError, "test failure"):
                thickness.main("input.json")
        gateway.dispose.assert_called_once_with()

    def test_orientation_and_summary_on_projection(self):
        with tempfile.TemporaryDirectory(prefix="uma projection ") as folder:
            root = Path(folder)
            images = root / "Images"
            analysis = root / "Analysis"
            images.mkdir()
            analysis.mkdir()
            y, x = np.mgrid[:64, :64]
            projection = (220 * np.exp(-((x - 25) / 3) ** 2) *
                          np.exp(-((y - 32) / 20) ** 2)).astype(np.uint8)
            tifffile.imwrite(root / "sample_processed.tif", projection)
            (root / "._sample_processed.tif").write_bytes(b"macOS metadata")
            alignment.process_part2_orientationpy(str(root), str(images))
            (root / "Tables" / "._ignored.csv").write_text("Not CSV")
            alignment.process_part3(str(root), str(analysis), 15,
                                    {"sample_processed": {"number_of_z_stacks": 7,
                                                          "z_stack_type": "slices"}})
            data = pd.read_csv(analysis / "Alignment_Summary.csv")
            self.assertEqual(len(data), 1)
            self.assertEqual(int(data.iloc[0]["Number_of_Z_Stacks"]), 7)
            value = data.iloc[0]["Percentage_Fibers_Aligned_Within_15_Degree"]
            self.assertTrue(np.isfinite(value))
            self.assertGreaterEqual(value, 0)
            self.assertLessEqual(value, 100)
            self.assertEqual(len(list(images.glob("*.png"))), 1)
            self.assertEqual(len(list((images / "normalized_images").glob("*.png"))), 1)


if __name__ == "__main__":
    unittest.main()
