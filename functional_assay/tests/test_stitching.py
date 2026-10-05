"""Stitching input, grid, output replacement, and command cleanup checks."""

import contextlib
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from functional_assay import stitching
from functional_assay.workflow import assay_directory
from uma_tools import cli


class StitchingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="uma stitch ")
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        calibration = patch.object(
            stitching,
            "read_well_calibration",
            return_value={
                "pixel_size_x_um": 0.5,
                "pixel_size_y_um": 0.25,
                "pixel_area_um2": 0.125,
            },
        )
        calibration.start()
        self.addCleanup(calibration.stop)

    def frames(self, well="WellA1", indices=range(9), prefix="sample"):
        files = []
        for index in indices:
            path = self.folder / (
                f"{prefix}__{well}_PointA1_{index:04d}_ChannelGFP.nd2"
            )
            path.write_text(str(index))
            files.append((index, path))
        return files

    def test_only_visible_image_files_are_grouped(self):
        first = self.frames()
        self.frames(well="WellB2")
        (self.folder / ("._" + first[0][1].name)).write_text("metadata")
        (self.folder / ("." + first[0][1].name)).write_text("hidden")
        (self.folder / "extra__WellA1_PointA1_0000_ChannelGFP.nd2").mkdir()
        (self.folder / "notes.txt").write_text("not an image")
        groups = stitching.discover_wells(self.folder)
        self.assertEqual(set(groups), {"WellA1", "WellB2"})
        self.assertEqual(
            [index for index, _ in groups["WellA1"]], list(range(9))
        )
        for files in groups.values():
            stitching.validate_frames(files)

    def test_nine_files_with_duplicate_and_missing_frame_are_rejected(self):
        files = self.frames(indices=range(8))
        files += self.frames(indices=[7], prefix="repeat")
        with self.assertRaisesRegex(
            ValueError, r"missing=\[8\].*duplicates=\[7\]"
        ):
            stitching.validate_frames(files)

    def test_missing_and_out_of_range_indices_are_rejected(self):
        for indices in (range(8), (*range(8), 9), range(10)):
            with self.subTest(indices=indices):
                with self.assertRaises(ValueError):
                    stitching.validate_frames(self.frames(indices=indices))

    def test_invalid_well_is_skipped_while_valid_well_is_stitched(self):
        self.frames()
        self.frames(well="WellB2", indices=range(8))
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            with patch.object(stitching, "stitch_well") as fuse:
                count = stitching.process_folder(self.folder, 32.8)
        self.assertEqual(count["completed_wells"], 1)
        self.assertEqual(count["status"], "PARTIAL")
        self.assertEqual(fuse.call_args.args[1], "WellA1")
        self.assertIn("Skipping WellB2", stream.getvalue())
        self.assertIn("missing=[8]", stream.getvalue())

    def test_previous_results_are_replaced_but_sources_are_retained(self):
        frames = self.frames()
        output = assay_directory(self.folder) / "Stitched_Results"
        (output / "old_nested").mkdir(parents=True)
        (output / "old_nested" / "old.tif").write_bytes(b"old")
        with patch.object(stitching, "stitch_well"):
            stitching.process_folder(self.folder, 32.8)
        self.assertFalse((output / "old_nested").exists())
        metadata = json.loads((output / "stitching_metadata.json").read_text())
        self.assertEqual(metadata["overlap_percent"], 32.8)
        self.assertEqual(metadata["status"], "SUCCESS")
        self.assertEqual(
            [
                item["grid_position"]
                for item in metadata["wells"]["WellA1"]["source_frames"]
            ],
            [1, 2, 3, 6, 5, 4, 7, 8, 9],
        )
        self.assertTrue(all(path.is_file() for _, path in frames))

    def test_invalid_attempt_removes_stale_results_and_records_failure(self):
        self.frames(indices=range(8))
        output = assay_directory(self.folder) / "Stitched_Results"
        output.mkdir()
        previous = output / "WellA1_stitched.tif"
        previous.write_bytes(b"old valid output")
        with patch.object(stitching, "stitch_well") as fuse:
            result = stitching.process_folder(self.folder, 32.8)
            self.assertEqual(result["completed_wells"], 0)
            self.assertEqual(result["status"], "FAILED")
        fuse.assert_not_called()
        self.assertFalse(previous.exists())

    def test_output_symlink_does_not_delete_its_target(self):
        target = self.folder / "originals"
        target.mkdir()
        original = target / "keep.nd2"
        original.write_bytes(b"keep")
        (assay_directory(self.folder) / "Stitched_Results").symlink_to(target)
        with self.assertRaisesRegex(ValueError, "symbolic link"):
            stitching.reset_output_folder(self.folder)
        self.assertEqual(original.read_bytes(), b"keep")

    def fake_imagej(self, channels=1, save_error=None):
        image = Mock()
        image.getNChannels.return_value = channels
        image.getWidth.return_value = 1200
        image.getHeight.return_value = 1200
        image.getNSlices.return_value = 3
        ij = Mock()
        events = []

        def run_macro(macro, input_image):
            self.assertIsNone(input_image)
            self.assertTrue(
                macro.startswith('run("Grid/Collection stitching", ')
            )
            options = json.loads(macro.split(", ", 1)[1][:-2])
            directory = Path(re.search(r"directory=\[(.*?)\]", options)[1])
            self.assertEqual(
                [
                    int((directory / f"image_{i}.nd2").read_text())
                    for i in range(1, 10)
                ],
                [0, 1, 2, 5, 4, 3, 6, 7, 8],
            )
            self.assertIn("fusion_method=[Linear Blending]", options)
            self.assertIn("tile_overlap=32.8", options)
            self.assertIn("grid_size_x=3 grid_size_y=3", options)
            self.assertIn("regression_threshold=0.30", options)
            self.assertIn("max/avg_displacement_threshold=2.50", options)
            self.assertIn("absolute_displacement_threshold=3.50", options)
            self.assertNotIn("compute_overlap", options)
            events.append("fused")
            return image

        def run_filter(received, name, options):
            self.assertIs(received, image)
            self.assertEqual((name, options), ("Sharpen", "stack"))
            events.append("sharpened")

        def save(received, kind, path):
            self.assertIs(received, image)
            self.assertEqual(kind, "Tiff")
            self.assertEqual(events, ["fused", "sharpened"])
            if save_error:
                raise save_error
            Path(path).write_bytes(b"test TIFF output")
            events.append("saved")

        interpreter = Mock()
        interpreter.return_value.runBatchMacro.side_effect = run_macro
        ij.run.side_effect = run_filter
        ij.saveAs.side_effect = save
        classes = {"ij.IJ": ij, "ij.macro.Interpreter": interpreter}
        return classes, image, events

    def test_metadata_records_real_dimensions_and_binds_output_pixels(self):
        from uma_tools.files import sha256_file

        files = self.frames()
        output = stitching.reset_output_folder(self.folder)
        classes, _, _ = self.fake_imagej()
        record = {}
        with patch("scyjava.jimport", side_effect=classes.__getitem__):
            path = stitching.stitch_well(
                self.folder, "WellA1", files, output, 32.8, record
            )
        self.assertEqual(record["width_px"], 1200)
        self.assertEqual(record["height_px"], 1200)
        self.assertEqual(record["slices"], 3)
        self.assertEqual(record["sha256"], sha256_file(path))

    def test_grid_sharpen_tiff_and_temporary_file_cleanup(self):
        files = self.frames()
        output = stitching.reset_output_folder(self.folder)
        classes, image, events = self.fake_imagej()
        with patch("scyjava.jimport", side_effect=classes.__getitem__):
            result = stitching.stitch_well(
                self.folder, "WellA1", files, output, 32.8
            )
        self.assertEqual(result.name, "WellA1_stitched.tif")
        self.assertEqual(events, ["fused", "sharpened", "saved"])
        self.assertEqual(list(output.iterdir()), [result])
        self.assertEqual(list(self.folder.glob(".uma_stitch_*")), [])
        image.close.assert_called_once_with()

    def test_copy_fallback_preserves_sources_and_cleans_up(self):
        files = self.frames()
        output = stitching.reset_output_folder(self.folder)
        classes, image, _ = self.fake_imagej()
        with patch("scyjava.jimport", side_effect=classes.__getitem__):
            with patch.object(
                Path, "symlink_to", side_effect=OSError("no symlinks")
            ):
                stitching.stitch_well(
                    self.folder, "WellA1", files, output, 32.8
                )
        self.assertTrue(all(path.exists() for _, path in files))
        self.assertEqual(list(self.folder.glob(".uma_stitch_*")), [])
        image.close.assert_called_once_with()

    def test_save_failure_closes_image_and_removes_temporary_files(self):
        files = self.frames()
        output = stitching.reset_output_folder(self.folder)
        classes, image, _ = self.fake_imagej(save_error=OSError("disk full"))
        with patch("scyjava.jimport", side_effect=classes.__getitem__):
            with self.assertRaisesRegex(OSError, "disk full"):
                stitching.stitch_well(
                    self.folder, "WellA1", files, output, 32.8
                )
        image.close.assert_called_once_with()
        self.assertEqual(list(self.folder.glob(".uma_stitch_*")), [])
        self.assertTrue(all(path.exists() for _, path in files))

    def test_multichannel_result_is_rejected_before_filter_or_save(self):
        files = self.frames()
        output = stitching.reset_output_folder(self.folder)
        classes, image, events = self.fake_imagej(channels=2)
        with patch("scyjava.jimport", side_effect=classes.__getitem__):
            with self.assertRaisesRegex(ValueError, "one channel"):
                stitching.stitch_well(
                    self.folder, "WellA1", files, output, 32.8
                )
        self.assertEqual(events, ["fused"])
        image.close.assert_called_once_with()
        self.assertEqual(list(output.iterdir()), [])

    def test_metadata_json_is_rejected_before_reading_or_starting_fiji(self):
        with patch.object(stitching, "initialize_imagej") as start:
            with patch.object(Path, "read_text") as read:
                with self.assertRaisesRegex(ValueError, "metadata"):
                    stitching.process_wells_stitching(
                        str(self.folder / "._input.json")
                    )
        start.assert_not_called()
        read.assert_not_called()

    def test_cli_overlap_and_worker_cleanup(self):
        for extra, expected in (([], 32.8), (["--overlap", "25"], 25.0)):
            with self.subTest(extra=extra):
                with patch.object(
                    stitching, "process_wells_stitching"
                ) as process:
                    with patch.object(
                        cli, "_shutdown_imagej_workers"
                    ) as shutdown:
                        self.assertEqual(
                            stitching.main(["-i", "input.json", *extra]), 0
                        )
                process.assert_called_once_with("input.json", expected)
                shutdown.assert_called_once_with()

    def test_invalid_overlap_is_rejected_before_analysis(self):
        for value in ("-1", "100", "nan", "inf", "abc"):
            with self.subTest(value=value):
                with patch.object(
                    stitching, "process_wells_stitching"
                ) as process:
                    with contextlib.redirect_stderr(io.StringIO()):
                        with self.assertRaises(SystemExit) as error:
                            stitching.main(
                                ["-i", "input.json", "--overlap", value]
                            )
                self.assertEqual(error.exception.code, 2)
                process.assert_not_called()

    def test_context_and_workers_close_on_success_and_processing_error(self):
        configuration = self.folder / "input.json"
        configuration.write_text(
            json.dumps({"folder_paths": [str(self.folder)]})
        )
        for error in (None, RuntimeError("stitching failed")):
            with self.subTest(error=error):
                context = Mock()

                def process(folder, overlap, ensure):
                    ensure()
                    if error:
                        raise error
                    return {"status": "SUCCESS"}

                with patch.object(
                    stitching, "initialize_imagej", return_value=context
                ):
                    with patch.object(
                        stitching,
                        "process_folder",
                        side_effect=process,
                    ):
                        with patch("scyjava.jvm_started", return_value=False):
                            with patch("scyjava.config.add_option") as option:
                                with patch.object(
                                    cli, "_shutdown_imagej_workers"
                                ) as workers:
                                    self.assertEqual(
                                        stitching.main(
                                            ["-i", str(configuration)]
                                        ),
                                        1 if error else 0,
                                    )
                option.assert_called_once_with("-Xmx16g")
                context.dispose.assert_called_once_with()
                workers.assert_called_once_with()

    def test_help_and_version_do_not_start_fiji(self):
        for option in ("--help", "--version"):
            script = (
                "import sys; from functional_assay.stitching import main\n"
                f"try: main([{option!r}])\n"
                "except SystemExit as error: assert error.code == 0\n"
                "assert 'imagej' not in sys.modules\n"
                "assert 'scyjava' not in sys.modules\n"
            )
            result = subprocess.run(
                [sys.executable, "-c", script],
                cwd=self.folder,
                capture_output=True,
                text=True,
                timeout=20,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

    @unittest.skipUnless(
        os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
        "Set UMA_RUN_IMAGEJ_TESTS=1 for the Java process-exit check",
    )
    def test_command_exits_after_java_work_on_success_and_failure(self):
        script = textwrap.dedent("""
            import os
            import sys
            import jpype
            from functional_assay import cell_analysis, stitching
            from uma_tools import cli

            def java_work(*_):
                context = None
                jar = os.environ.get("UMA_TEST_IMAGEJ_JAR")
                if jar:
                    jpype.startJVM("-Djava.awt.headless=true", classpath=[jar])
                else:
                    from uma_tools.imagej import initialize_imagej
                    context = initialize_imagej()
                try:
                    threads = jpype.JClass("ij.util.ThreadUtil")
                    pool = threads.threadPoolExecutor
                    task = jpype.JProxy(
                        "java.lang.Runnable", dict(run=lambda: None)
                    )
                    pool.submit(task).get()
                    assert pool.getPoolSize() > 0
                    print("Java worker was started", flush=True)
                    if os.environ["UMA_STITCH_TEST_FAILS"] == "1":
                        raise RuntimeError("deliberate stitching failure")
                finally:
                    if context is not None:
                        context.dispose()

            if os.environ["UMA_TEST_COMMAND"] == "uma_cell_count":
                cell_analysis.run_analysis = java_work
                sys.exit(cell_analysis.main(["-i", "input.json"]))
            else:
                stitching.process_wells_stitching = java_work
                sys.exit(stitching.main(["-i", "input.json"]))
        """)
        for command, fail in (
            ("uma_stitching", False),
            ("uma_stitching", True),
            ("uma_cell_count", False),
            ("uma_cell_count", True),
        ):
            with self.subTest(command=command, analysis_fails=fail):
                result = subprocess.run(
                    [sys.executable, "-c", script],
                    env=dict(
                        os.environ,
                        UMA_STITCH_TEST_FAILS=str(int(fail)),
                        UMA_TEST_COMMAND=command,
                    ),
                    cwd=self.folder,
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                self.assertIn("Java worker was started", result.stdout)
                self.assertEqual(
                    result.returncode, 1 if fail else 0, result.stderr
                )
                if fail:
                    self.assertIn(
                        "deliberate stitching failure", result.stderr
                    )


if __name__ == "__main__":
    unittest.main()
