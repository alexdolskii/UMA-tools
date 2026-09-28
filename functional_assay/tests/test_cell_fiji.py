"""Opt-in comparisons against native Fiji macros and ParticleAnalyzer."""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from functional_assay import cell_imagej


@unittest.skipUnless(
    os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
    "Set UMA_RUN_IMAGEJ_TESTS=1 for native Fiji validation",
)
class NativeFijiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The optional archive bootstrap is for offline tests only. Production
        # always uses UMA's shared, pinned Fiji endpoint.
        import imagej
        import jpype
        from scyjava import jimport

        from uma_tools.imagej import initialize_imagej

        assert imagej is not None
        directory = os.environ.get("UMA_TEST_FIJI_DIR")
        if directory:
            fiji = Path(directory)
            jpype.startJVM(
                "-Xmx4g",
                "-Djava.awt.headless=true",
                f"-Dplugins.dir={fiji / 'plugins'}",
                "--add-opens=java.base/java.lang=ALL-UNNAMED",
                "--add-opens=java.base/java.util=ALL-UNNAMED",
                "--add-opens=java.desktop/sun.awt.X11=ALL-UNNAMED",
                classpath=[str(path) for path in fiji.rglob("*.jar")],
            )
            jpype.JClass("net.imagej.patcher.LegacyInjector").preinit()
        cls.context = initialize_imagej()
        cls.ij = jimport("ij.IJ")
        cls.jimport = staticmethod(jimport)
        jimport("ij.Prefs").blackBackground = True

    @classmethod
    def tearDownClass(cls):
        from uma_tools.imagej import shutdown_imagej_workers

        try:
            cls.context.dispose()
        finally:
            shutdown_imagej_workers()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="UMA native ")
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def native_particles(self, mask, minimum, exclude=False):
        analyzer_class = self.jimport("ij.plugin.filter.ParticleAnalyzer")
        table = self.jimport("ij.measure.ResultsTable")()
        options = analyzer_class.SHOW_MASKS
        if exclude:
            options |= analyzer_class.EXCLUDE_EDGE_PARTICLES
        analyzer = analyzer_class(options, 1, table, minimum, float("inf"))
        analyzer.setHideOutputImage(True)
        image = cell_imagej.mask_image(mask, "Reference particles")
        image.getProcessor().setThreshold(255, 255, 0)
        output = None
        try:
            self.assertTrue(analyzer.analyze(image))
            output = analyzer.getOutputImage()
            pixels = np.asarray(output.getProcessor().getPixels())
            accepted = pixels.reshape(mask.shape) != 0
            areas = [table.getValue("Area", i) for i in range(table.size())]
            return accepted, np.array(areas)
        finally:
            image.close()
            if output is not None:
                output.close()

    def test_pixel_connectivity_holes_sizes_and_edges_match_imagej(self):
        mask = np.zeros((32, 32), dtype=bool)
        mask[6:12, 6:12] = True
        mask[7:11, 7:11] = False
        mask[0:2, 0:3] = True
        mask[20, 20] = mask[21, 21] = True
        mask[26, 26] = True
        for minimum in (2, 5, 5.1):
            for exclude in (False, True):
                with self.subTest(minimum=minimum, exclude=exclude):
                    native_mask, native_areas = self.native_particles(
                        mask, minimum, exclude
                    )
                    labels, areas = cell_imagej.filter_particles(
                        mask, minimum, exclude_edges=exclude
                    )
                    np.testing.assert_array_equal(labels != 0, native_mask)
                    np.testing.assert_array_equal(
                        np.sort(areas), np.sort(native_areas)
                    )

    def fixture(self):
        import tifffile

        rng = np.random.default_rng(384)
        stack = rng.integers(10, 20, size=(3, 128, 128), dtype=np.uint16)
        y, x = np.ogrid[:128, :128]
        for z, cy, cx, radius, intensity in (
            (0, 35, 35, 12, 3000),
            (1, 35, 53, 12, 2600),
            (2, 85, 85, 9, 4000),
            (0, 3, 100, 12, 3200),
            (1, 100, 20, 1, 1000),
        ):
            stack[z][(y - cy) ** 2 + (x - cx) ** 2 <= radius**2] = intensity
        path = self.folder / "WellA1_stitched.tif"
        tifffile.imwrite(path, stack, imagej=True, metadata={"axes": "ZYX"})
        return path

    def reference(self, path, thresholds):
        source = self.ij.openImage(str(path))
        macro = (
            'run("Z Project...", "projection=[Max Intensity]");'
            'run("Subtract Background...", "rolling=50");'
            'run("Enhance Contrast", "saturated=0.35");'
            'run("Smooth"); run("Smooth");'
        )
        projection = self.jimport("ij.macro.Interpreter")().runBatchMacro(
            macro, source
        )
        try:
            if thresholds is None:
                self.ij.setAutoThreshold(
                    projection, "RenyiEntropy dark no-reset"
                )
            else:
                self.ij.setThreshold(projection, *thresholds)
            processor = projection.getProcessor()
            shape = (projection.getHeight(), projection.getWidth())
            pixels = (
                np.asarray(processor.getPixels())
                .astype(np.uint16)
                .reshape(shape)
            )
            raw_mask = (
                np.asarray(processor.createMask().getPixels()).reshape(shape)
                != 0
            )
            area_mask, _ = self.native_particles(raw_mask, 5)
            divided = cell_imagej.mask_image(area_mask, "Native watershed")
            try:
                self.ij.run(divided, "Watershed", "")
                watershed = (
                    np.asarray(divided.getProcessor().getPixels()).reshape(
                        shape
                    )
                    != 0
                )
                counted, areas = self.native_particles(watershed, 5, True)
            finally:
                divided.close()
            return (
                pixels,
                area_mask,
                counted,
                areas,
                (
                    float(processor.getMinThreshold()),
                    float(processor.getMaxThreshold()),
                ),
            )
        finally:
            projection.close()
            source.close()

    def test_processing_and_masks_match_native_reference(self):
        path = self.fixture()
        for thresholds in (None, (50, 65535), (100, 2500)):
            with self.subTest(thresholds=thresholds):
                pixels, area, count, sizes, bounds = self.reference(
                    path, thresholds
                )
                result = cell_imagej.analyze_image(path, thresholds, 5)
                np.testing.assert_array_equal(result["projection"], pixels)
                np.testing.assert_array_equal(result["area_mask"], area)
                np.testing.assert_array_equal(result["counting_mask"], count)
                np.testing.assert_array_equal(
                    np.sort(result["object_areas_px2"]), np.sort(sizes)
                )
                self.assertEqual(
                    (result["threshold_lower"], result["threshold_upper"]),
                    bounds,
                )
                self.assertGreater(int(area.sum()), int(count.sum()))

    def test_calibrated_masks_round_trip_and_qc_has_true_contours(self):
        from PIL import Image

        result = cell_imagej.analyze_image(self.fixture(), None, 5)
        path = self.folder / "mask.tif"
        cell_imagej.save_mask(
            path,
            result["area_mask"],
            {
                "pixel_size_x_um": 0.5,
                "pixel_size_y_um": 0.25,
            },
        )
        saved = self.ij.openImage(str(path))
        try:
            scale = saved.getCalibration()
            self.assertAlmostEqual(scale.pixelWidth, 0.5)
            self.assertAlmostEqual(scale.pixelHeight, 0.25)
            self.assertIn(str(scale.getUnit()), ("µm", "um", "micron"))
            recovered = (
                np.asarray(saved.getProcessor().getPixels()).reshape(
                    result["area_mask"].shape
                )
                != 0
            )
            np.testing.assert_array_equal(recovered, result["area_mask"])
        finally:
            saved.close()
        cell_imagej.save_contours(self.folder / "qc.png", result)
        with Image.open(self.folder / "qc.png") as image:
            rgb = np.asarray(image)
            green = np.all(rgb == (0, 255, 0), axis=2)
        self.assertTrue(green.any())
        self.assertTrue(np.all(result["counting_mask"][green]))

    def test_blank_image_and_invalid_depth(self):
        import tifffile

        path = self.folder / "blank_stitched.tif"
        tifffile.imwrite(path, np.zeros((32, 32), dtype=np.uint16))
        result = cell_imagej.analyze_image(path, (50, 65535), 5)
        self.assertEqual(result["area_mask"].sum(), 0)
        self.assertEqual(len(result["object_areas_px2"]), 0)
        tifffile.imwrite(path, np.zeros((32, 32), dtype=np.uint8))
        with self.assertRaisesRegex(ValueError, "16-bit"):
            cell_imagej.analyze_image(path, None, 5)


if __name__ == "__main__":
    unittest.main()
