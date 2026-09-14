"""Opt-in Fiji integration checks on synthetic TIFF stacks, without a display.

Run with UMA_RUN_IMAGEJ_TESTS=1 after installing the Conda environment and package.
These checks do not replace validation on representative experimental ND2 files.
"""

import json
import logging
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock


@unittest.skipUnless(os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
                     "Set UMA_RUN_IMAGEJ_TESTS=1 to download and run Fiji")
class HeadlessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import imagej
        cls.ij = imagej.init("sc.fiji:fiji:2.14.0", mode="headless")

    @classmethod
    def tearDownClass(cls):
        from scyjava import jimport
        try:
            cls.ij.dispose()
        finally:
            # ImageJ1 filters leave non-daemon workers outside the Fiji context.
            # Close the shared pool only after the final integration test.
            pool = jimport("ij.util.ThreadUtil").threadPoolExecutor
            pool.shutdown()
            seconds = jimport("java.util.concurrent.TimeUnit").SECONDS
            if not pool.awaitTermination(5, seconds):
                raise RuntimeError("ImageJ worker pool did not terminate")

    def setUp(self):
        import numpy as np
        import tifffile
        self.temp = tempfile.TemporaryDirectory(prefix="uma assay test ")
        self.folder = Path(self.temp.name)
        stack = np.zeros((7, 32, 32), dtype=np.uint16)
        stack[1:6, 5:27, 8:12] = 30000
        stack[2:5, 9:23, 21:25] = 45000
        tifffile.imwrite(self.folder / "sample.tiff", stack, imagej=True,
                         resolution=(2, 2),
                         metadata={"axes": "ZYX", "spacing": 1.0, "unit": "um"})
        (self.folder / "._sample.tiff").write_bytes(b"Not an image: macOS metadata")
        self.manifest = self.folder / "input.json"
        self.manifest.write_text(json.dumps({"folder_paths": [str(self.folder)]}))
        self.handlers = list(logging.getLogger().handlers)

    def tearDown(self):
        for handler in list(logging.getLogger().handlers):
            if handler not in self.handlers:
                logging.getLogger().removeHandler(handler)
                handler.close()
        self.temp.cleanup()

    def test_alignment_outputs_and_metadata_filter(self):
        import pandas as pd
        from uma_tools import alignment_analysis as alignment
        self.assertEqual(alignment.get_folder_paths(str(self.manifest)), [str(self.folder)])
        alignment.process_folder(str(self.folder), 1, 15, 64, 64, self.ij)
        summaries = list(self.folder.glob("Alignment_assay_results_*/Analysis/Alignment_Summary.csv"))
        self.assertEqual(len(summaries), 1)
        data = pd.read_csv(summaries[0])
        self.assertEqual(len(data), 1)
        self.assertIn("sample", data.iloc[0]["File_Name"])
        self.assertEqual(int(data.iloc[0]["Number_of_Z_Stacks"]), 7)
        self.assertTrue(list(self.folder.glob("Alignment_assay_results_*/Images/*.png")))

    def test_thickness_outputs_and_metadata_filter(self):
        import numpy as np
        import pandas as pd
        from uma_tools import thickness_analysis as thickness
        ij, image_plus, windows, table, duplicator, system = thickness.import_java_classes()
        thickness.process_single_folder(self.ij, ij, windows, duplicator, table,
                                        str(self.folder), ".tiff", 1)
        summaries = list(self.folder.glob("Thickness_assay_results_*/Thickness_Summary.csv"))
        self.assertEqual(len(summaries), 1)
        data = pd.read_csv(summaries[0])
        self.assertEqual(len(data), 1)
        self.assertEqual(list(data.columns), ["File_Name", "Area", "StdDev", "Min", "Max", "Median"])
        self.assertTrue(np.isfinite(data.iloc[0, 1:].astype(float)).all())
        self.assertGreater(float(data.iloc[0]["Median"]), 0)
        self.assertEqual(len(list(summaries[0].parent.glob("Mask_*.tif"))), 1)
        self.assertEqual(len(list(summaries[0].parent.glob("Local_Thickness_*.tif"))), 1)

    def test_local_thickness_api_matches_original_menu_command(self):
        import numpy as np
        from scyjava import jimport
        ij = jimport("ij.IJ")
        image = ij.createImage("mask", "8-bit black", 24, 24, 1)
        image.getProcessor().setValue(255)
        image.getProcessor().setRoi(6, 6, 12, 12)
        image.getProcessor().fill()
        image.getProcessor().resetRoi()
        image.getCalibration().pixelWidth = 0.5
        image.getCalibration().pixelHeight = 0.5
        direct_input = image.duplicate()
        reference = jimport("ij.macro.Interpreter")().runBatchMacro(
            'run("Local Thickness (masked, calibrated, silent)");', image)
        plugin = jimport("sc.fiji.localThickness.LocalThicknessWrapper")()
        plugin.setSilence(True)
        plugin.setShowOptions(False)
        plugin.maskThicknessMap = True
        plugin.calibratePixels = True
        actual = plugin.processImage(direct_input)
        np.testing.assert_allclose(np.array(actual.getProcessor().getPixels()),
                                   np.array(reference.getProcessor().getPixels()),
                                   rtol=0, atol=0, equal_nan=True)
        self.assertEqual(actual.getCalibration().pixelWidth,
                         reference.getCalibration().pixelWidth)
        for imp in (actual, reference, direct_input, image):
            imp.close()

    def test_thickness_matches_original_projection_with_calibration(self):
        import numpy as np
        import tifffile
        from uma_tools import thickness_analysis as thickness

        filename = "calibrated.tiff"
        stack = np.zeros((17, 64, 64), dtype=np.uint16)
        stack[3:14, 8:56, 10:28] = 30000
        stack[5:12, 12:52, 37:55] = 45000
        tifffile.imwrite(self.folder / filename, stack, imagej=True,
                         resolution=(2, 2),
                         metadata={"axes": "ZYX", "spacing": 0.75, "unit": "um"})

        IJ, image_plus, windows, table, duplicator, _ = thickness.import_java_classes()
        loaded = self.ij.convert().convert(
            self.ij.io().open(str(self.folder / filename)), image_plus)
        try:
            self.assertAlmostEqual(loaded.getCalibration().pixelWidth, 0.5)
            self.assertAlmostEqual(loaded.getCalibration().pixelDepth, 0.75)
        finally:
            loaded.close()

        real_jimport = thickness.jimport

        class OriginalMacroProjector:
            """Reference adapter for both the former and current Python call sites."""
            MAX_METHOD = 1

            def __init__(self, image):
                self.image = image
                self.projection = None

            def setMethod(self, method):
                if method != self.MAX_METHOD:
                    raise AssertionError("The reference requires maximum projection")

            def doProjection(self, *args):
                self.projection = self.run(self.image, "max")

            def getProjection(self):
                return self.projection

            @staticmethod
            def run(image, method):
                if method != "max":
                    raise AssertionError("The reference requires maximum projection")
                return real_jimport("ij.macro.Interpreter")().runBatchMacro(
                    'run("Z Project...", "projection=[Max Intensity]");', image)

        actual_dir = self.folder / "actual"
        reference_dir = self.folder / "reference"
        actual_dir.mkdir()
        reference_dir.mkdir()
        actual = thickness.process_single_file(
            self.ij, IJ, windows, duplicator, table, str(self.folder),
            filename, 1, str(actual_dir))
        with mock.patch.object(
                thickness, "jimport", side_effect=lambda name:
                OriginalMacroProjector if name == "ij.plugin.ZProjector"
                else real_jimport(name)):
            reference = thickness.process_single_file(
                self.ij, IJ, windows, duplicator, table, str(self.folder),
                filename, 1, str(reference_dir))

        self.assertIsNotNone(actual)
        self.assertIsNotNone(reference)
        for prefix in ("Mask_", "Local_Thickness_"):
            actual_image = IJ.openImage(str(actual_dir / f"{prefix}{filename}.tif"))
            reference_image = IJ.openImage(str(reference_dir / f"{prefix}{filename}.tif"))
            try:
                self.assertIsNotNone(actual_image)
                self.assertIsNotNone(reference_image)
                self.assertEqual(tuple(actual_image.getDimensions()),
                                 tuple(reference_image.getDimensions()))
                actual_cal = actual_image.getCalibration()
                reference_cal = reference_image.getCalibration()
                self.assertAlmostEqual(reference_cal.pixelWidth, 0.75)
                self.assertAlmostEqual(reference_cal.pixelHeight, 0.5)
                self.assertEqual(actual_cal.getUnit(), reference_cal.getUnit())
                for axis in ("pixelWidth", "pixelHeight", "pixelDepth"):
                    self.assertEqual(getattr(actual_cal, axis),
                                     getattr(reference_cal, axis))
                np.testing.assert_array_equal(
                    np.array(actual_image.getProcessor().getPixels()),
                    np.array(reference_image.getProcessor().getPixels()))
            finally:
                for image in (actual_image, reference_image):
                    if image is not None:
                        image.close()

        self.assertEqual(list(actual), list(reference))
        self.assertEqual(actual["File_Name"], reference["File_Name"])
        columns = ("Area", "StdDev", "Min", "Max", "Median")
        self.assertTrue(np.isfinite([reference[name] for name in columns]).all())
        self.assertGreater(reference["Max"], 0)
        np.testing.assert_allclose([actual[name] for name in columns],
                                   [reference[name] for name in columns],
                                   rtol=0, atol=0)


if __name__ == "__main__":
    unittest.main()
