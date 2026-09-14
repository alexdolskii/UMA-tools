"""Behavioral contracts for modular assays and their command entry points."""

import contextlib
import csv
import hashlib
import io
import json
import logging
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import tifffile

from uma_tools import collection as collector
from uma_tools.assays import alignment, thickness

FIXTURES = Path(__file__).with_name("fixtures")
REPOSITORY = Path(__file__).resolve().parents[1]
ENTRY_POINTS = (
    ("1_alignment.py", "uma_alignment"),
    ("2_thickness.py", "uma_thickness"),
    ("3_area.py", "area_analysis"),
    ("4_collect_results.py", "uma_collect_results"),
    ("5_report.py", "uma_report"),
)


def crossing_projection():
    """Return the deterministic uint8 image used by the 0.2.5 reference."""
    y, x = np.mgrid[:64, :64]
    first = 210 * np.exp(-(((x - 0.36 * y - 15) / 2.8) ** 2))
    first *= np.exp(-(((y - 32) / 24) ** 2))
    second = 100 * np.exp(-(((x + 0.8 * y - 74) / 3.5) ** 2))
    second *= np.exp(-(((y - 32) / 25) ** 2))
    return np.minimum(255, first + second).astype(np.uint8)


class FrozenAlignmentTests(unittest.TestCase):
    def test_distribution_summary_and_rgb_match_frozen_025_baseline(self):
        """Compare real calculations with the unmodified release outputs."""
        expected = json.loads(
            (FIXTURES / "alignment_025_crossing.json").read_text(
                encoding="utf-8"
            )
        )
        projection = crossing_projection()
        self.assertEqual(
            hashlib.sha256(projection.tobytes()).hexdigest(),
            expected["projection_sha256"],
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            images = root / "Images"
            analysis = root / "Analysis"
            images.mkdir()
            analysis.mkdir()
            tifffile.imwrite(root / "crossing_processed.tif", projection)
            sampled_rgb = []
            coordinates = np.array(expected["rgb_sample_coordinates"])
            converter = alignment.matplotlib.colors.hsv_to_rgb

            def capture(values):
                rgb = converter(values)
                if values.shape == (64, 64, 3):
                    sampled_rgb.append(
                        rgb[coordinates[:, 0], coordinates[:, 1]].copy()
                    )
                return rgb

            with (
                contextlib.redirect_stdout(io.StringIO()),
                patch.object(
                    alignment.matplotlib.colors,
                    "hsv_to_rgb",
                    side_effect=capture,
                ),
            ):
                alignment.process_part2_orientationpy(str(root), str(images))
                alignment.process_part3(
                    str(root),
                    str(analysis),
                    15,
                    {
                        "crossing_processed": {
                            "number_of_z_stacks": 7,
                            "z_stack_type": "slices",
                        }
                    },
                )
            distribution = pd.read_csv(
                root
                / "Tables"
                / "crossing_processed_orientation_distribution.csv"
            )
            processed = pd.read_csv(
                analysis
                / "crossing_processed_orientation_distribution_processed.csv"
            )
            for actual, key in (
                (distribution, "distribution"),
                (processed, "processed_distribution"),
            ):
                self.assertEqual(list(actual), list(expected[key]))
                reference = pd.DataFrame(expected[key])
                for column in actual:
                    observed_values = actual[column].to_numpy()
                    reference_values = reference[column].to_numpy()
                    if column == "perc_occvalue2sum_of_occvalue":
                        # Division and decimal CSV parsing may differ by
                        # a few float64 ULPs between x86_64 and arm64.
                        # Counts, bins, ranks, and summary remain exact.
                        if not np.array_equal(
                            observed_values, reference_values
                        ):
                            difference = np.max(
                                np.abs(observed_values - reference_values)
                            )
                            print(
                                "Percentage CSV rounding: maximum absolute "
                                f"difference {difference:.17g}"
                            )
                        np.testing.assert_array_max_ulp(
                            observed_values, reference_values, maxulp=4
                        )
                    else:
                        np.testing.assert_array_equal(
                            observed_values, reference_values
                        )
            summary = pd.read_csv(analysis / "Alignment_Summary.csv")
            self.assertEqual(
                summary.to_dict(orient="records"), [expected["summary"]]
            )
            self.assertGreaterEqual(len(sampled_rgb), 2)
            # RGB arrays are float32; allow one small cross-platform rounding
            # difference. Histogram and angle columns stay exact.
            np.testing.assert_allclose(
                [sampled_rgb[0], sampled_rgb[-1]],
                [expected["rgb_samples"][0], expected["rgb_samples"][-1]],
                rtol=0,
                atol=1e-6,
            )
            self.assertEqual(len(list(images.glob("*.png"))), 1)
            self.assertEqual(
                len(list((images / "normalized_images").glob("*.png"))), 1
            )


class EntryPointTests(unittest.TestCase):
    def test_numbered_scripts_help_and_version_are_lightweight(
        self,
    ):
        driver = """
import runpy
import sys
path, argument = sys.argv[1:]
sys.argv = [path, argument]
try:
    runpy.run_path(path, run_name="__main__")
except SystemExit as error:
    assert error.code in (None, 0), error.code
for name in ("imagej", "jpype", "matplotlib.pyplot"):
    assert name not in sys.modules, name
"""
        with tempfile.TemporaryDirectory() as cwd:
            for filename, _ in ENTRY_POINTS:
                for argument in ("--help", "--version"):
                    with self.subTest(file=filename, argument=argument):
                        path = REPOSITORY / "code" / filename
                        completed = subprocess.run(
                            [
                                sys.executable,
                                "-c",
                                driver,
                                str(path),
                                argument,
                            ],
                            cwd=cwd,
                            capture_output=True,
                            text=True,
                            timeout=30,
                        )
                        self.assertEqual(
                            completed.returncode, 0, completed.stderr
                        )
                        self.assertTrue(completed.stdout.strip())

    def test_installed_commands_keep_argument_error_exit_status(self):
        with tempfile.TemporaryDirectory() as cwd:
            for _, command in ENTRY_POINTS:
                executable = Path(sys.executable).parent / command
                with self.subTest(command=command):
                    completed = subprocess.run(
                        [str(executable)],
                        cwd=cwd,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    self.assertEqual(completed.returncode, 2, completed.stderr)
                    self.assertIn("--input", completed.stderr)

    def test_numbered_commands_reach_canonical_analyses_without_adapters(self):
        """Exercise launch, argument routing, and exit status in isolation."""
        driver = """
import runpy
import sys
import types

path, module_name, function_name, result = sys.argv[1:]
for removed in (
    "alignment_analysis", "thickness_analysis", "area_analysis",
    "collect_results", "report", "report_rendering",
    "reporting.engine", "reporting.collection",
):
    sys.modules["uma_tools." + removed] = None
module = types.ModuleType(module_name)
calls = []

def calculate(*args):
    calls.append(args)
    assert sys.argv[1:] == ["-i", "routing path.json"]
    return int(result)

setattr(module, function_name, calculate)
sys.modules[module_name] = module
sys.argv = [path, "-i", "routing path.json"]
try:
    runpy.run_path(path, run_name="__main__")
except SystemExit as error:
    assert error.code == int(result), error.code
else:
    raise AssertionError("The launcher did not propagate its exit status")
expected = {
    "uma_tools.assays.alignment": ("routing path.json", 15),
    "uma_tools.assays.thickness": ("routing path.json",),
}.get(module_name, ())
assert calls == [expected], calls
for heavy in ("imagej", "matplotlib.pyplot"):
    assert heavy not in sys.modules, heavy
"""
        targets = (
            ("assays.alignment", "main_fibronectin_processing", 0),
            ("assays.thickness", "main", 0),
            ("assays.area", "main", 3),
            ("collection", "main", 3),
            ("reporting.workflow", "main", 3),
        )
        with tempfile.TemporaryDirectory() as cwd:
            for (filename, _), (module, entry, status) in zip(
                ENTRY_POINTS, targets
            ):
                with self.subTest(command=filename):
                    completed = subprocess.run(
                        [
                            sys.executable,
                            "-c",
                            driver,
                            str(REPOSITORY / "code" / filename),
                            f"uma_tools.{module}",
                            entry,
                            str(status),
                        ],
                        cwd=cwd,
                        capture_output=True,
                        text=True,
                        timeout=30,
                    )
                    self.assertEqual(completed.returncode, 0, completed.stderr)


class ScopedAssayTests(unittest.TestCase):
    def test_image_extension_directories_do_not_reach_image_openers(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "sample.nd2").touch()
            (root / "._sample.nd2").touch()
            (root / ".hidden.nd2").touch()
            (root / "directory.nd2").mkdir()
            output = root / "output"
            output.mkdir()
            ij = Mock()
            ij.openImage.return_value = None
            with (
                contextlib.redirect_stdout(io.StringIO()),
                patch.object(alignment.sj, "jimport", return_value=ij),
            ):
                alignment.process_part1(
                    str(root), str(output), 1, 64, 64, Mock()
                )
            ij.openImage.assert_called_once_with(str(root / "sample.nd2"))
            with (
                contextlib.redirect_stdout(io.StringIO()),
                patch.object(
                    thickness, "process_single_file", return_value=None
                ) as process,
            ):
                thickness.process_single_folder(
                    Mock(),
                    Mock(),
                    Mock(),
                    Mock(),
                    Mock(),
                    str(root),
                    ".nd2",
                    1,
                )
            self.assertEqual(process.call_count, 1)
            self.assertEqual(process.call_args.args[6], "sample.nd2")

    def test_alignment_handlers_are_closed_after_success_and_failure(self):
        self.check_scoped_log(alignment)

    def test_thickness_handlers_are_closed_after_success_and_failure(self):
        self.check_scoped_log(thickness)

    def test_same_instant_creates_distinct_collector_recognized_runs(self):
        frozen_time = datetime(2026, 9, 14, 12, 34, 56, 123456)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for module, name in (
                (alignment, "Alignment"),
                (thickness, "Thickness"),
            ):
                with self.subTest(assay=name):
                    source = root / name
                    source.mkdir()
                    (source / "sample.nd2").touch()

                    def save_alignment(results, analysis, angle, metadata):
                        target = Path(analysis) / "Alignment_Summary.csv"
                        with target.open(
                            "w", newline="", encoding="utf-8"
                        ) as stream:
                            writer = csv.writer(stream)
                            writer.writerow(
                                [
                                    "File_Name",
                                    "Number_of_Z_Stacks",
                                    "Z_Stack_Type",
                                    "Percentage_Fibers_Aligned_Within_15_Degree",
                                    "Orientation_Mode",
                                ]
                            )
                            writer.writerow(
                                [
                                    "sample_processed_orientation_distribution.csv",
                                    7,
                                    "slices",
                                    75.0,
                                    "aligned",
                                ]
                            )

                    with contextlib.ExitStack() as stack:
                        stack.enter_context(
                            contextlib.redirect_stdout(io.StringIO())
                        )
                        clock = stack.enter_context(
                            patch.object(module, "datetime")
                        )
                        clock.now.return_value = frozen_time
                        if module is alignment:
                            stack.enter_context(
                                patch.object(
                                    module, "process_part1", return_value={}
                                )
                            )
                            stack.enter_context(
                                patch.object(
                                    module, "process_part2_orientationpy"
                                )
                            )
                            stack.enter_context(
                                patch.object(
                                    module,
                                    "process_part3",
                                    side_effect=save_alignment,
                                )
                            )

                            def invoke():
                                return module.process_folder(
                                    str(source), 1, 15, 64, 64, Mock()
                                )
                        else:
                            stack.enter_context(
                                patch.object(
                                    module,
                                    "process_single_file",
                                    return_value={
                                        "File_Name": "sample.nd2",
                                        "Area": 12.5,
                                        "StdDev": 0.5,
                                        "Min": 1.0,
                                        "Max": 4.0,
                                        "Median": 2.0,
                                    },
                                )
                            )

                            def invoke():
                                return module.process_single_folder(
                                    Mock(),
                                    Mock(),
                                    Mock(),
                                    Mock(),
                                    Mock(),
                                    str(source),
                                    ".nd2",
                                    1,
                                )

                        invoke()
                        invoke()
                    directories = list(source.glob("*assay_results_*"))
                    self.assertEqual(len(directories), 2)
                    self.assertNotEqual(
                        directories[0].name, directories[1].name
                    )
                    assay = next(a for a in collector.ASSAYS if a.name == name)
                    records = []
                    # The independent collector recognizes both new names.
                    # Equal scientific timestamps remain explicitly ambiguous.
                    with self.assertRaisesRegex(
                        collector.ValidationError, "same latest timestamp"
                    ):
                        collector.select_latest(source, assay, records, Mock())
                    self.assertEqual(len(records), 2)
                    self.assertEqual(
                        {row["Status"] for row in records}, {"AMBIGUOUS"}
                    )

    def check_scoped_log(self, module):
        root_logger = logging.getLogger()
        original_root_handlers = list(root_logger.handlers)
        logger = module._LOGGER
        original_handlers = list(logger.handlers)
        original_level = logger.level
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            seen_handlers = []
            first_logs = []
            for fails in (False, True):
                source = root / (
                    "failed folder" if fails else "successful folder"
                )
                source.mkdir()
                (source / "sample.nd2").touch()
                message = "failure-only message" if fails else "first message"

                def calculate(*args, **kwargs):
                    seen_handlers.extend(
                        h
                        for h in logger.handlers
                        if isinstance(h, logging.FileHandler)
                    )
                    logger.info(message)
                    if fails:
                        raise RuntimeError("deliberate analysis failure")
                    return {}

                with contextlib.ExitStack() as stack:
                    stack.enter_context(
                        contextlib.redirect_stdout(io.StringIO())
                    )
                    if module is alignment:
                        stack.enter_context(
                            patch.object(
                                module, "process_part1", side_effect=calculate
                            )
                        )
                        stack.enter_context(
                            patch.object(module, "process_part2_orientationpy")
                        )
                        stack.enter_context(
                            patch.object(module, "process_part3")
                        )

                        def invoke():
                            return module.process_folder(
                                str(source), 1, 15, 64, 64, Mock()
                            )
                    else:
                        stack.enter_context(
                            patch.object(
                                module,
                                "process_single_file",
                                side_effect=calculate,
                            )
                        )

                        def invoke():
                            return module.process_single_folder(
                                Mock(),
                                Mock(),
                                Mock(),
                                Mock(),
                                Mock(),
                                str(source),
                                ".nd2",
                                1,
                            )

                    if fails:
                        with self.assertRaisesRegex(
                            RuntimeError, "deliberate analysis failure"
                        ):
                            invoke()
                    else:
                        invoke()
                self.assertEqual(root_logger.handlers, original_root_handlers)
                self.assertEqual(logger.handlers, original_handlers)
                self.assertEqual(logger.level, original_level)
                self.assertTrue(seen_handlers)
                self.assertTrue(all(h.stream is None for h in seen_handlers))
                logs = list(source.glob("*assay_results_*/*.log"))
                self.assertEqual(len(logs), 1)
                self.assertIn(message, logs[0].read_text(encoding="utf-8"))
                if not fails:
                    first_logs = logs
            for log in first_logs:
                self.assertNotIn(
                    "failure-only message", log.read_text(encoding="utf-8")
                )


if __name__ == "__main__":
    unittest.main()
