"""Native command recovery and process exit, isolated from the test JVM."""

import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


@unittest.skipUnless(
    os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1", "Requires the Fiji runtime"
)
class PartialFijiTests(unittest.TestCase):
    def test_python_and_real_java_use_the_same_managed_temp_directory(self):
        from uma_tools.runtime import run_command

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            runtime_home = root / "runtime with spaces"
            worker = root / "java_temp_check.py"
            evidence = root / "evidence.json"
            worker.write_text("""
import os, sys, tempfile
from pathlib import Path
from uma_tools.imagej import initialize_imagej, shutdown_imagej_workers
from uma_tools.runtime import record_completion, write_json
from scyjava import jimport
ij = initialize_imagej()
try:
    system = jimport("java.lang.System")
    file = jimport("java.io.File").createTempFile("uma-test-", ".tmp")
    write_json(Path(sys.argv[2]), {
        "python_temp": tempfile.gettempdir(),
        "java_temp": str(system.getProperty("java.io.tmpdir")),
        "java_file": str(file.getAbsolutePath()),
        "managed_temp": os.environ["TMPDIR"],
    })
finally:
    try:
        ij.dispose()
    finally:
        shutdown_imagej_workers()
record_completion(0, True)
""")
            with mock.patch.dict(
                os.environ, UMA_RUNTIME_HOME=str(runtime_home)
            ):
                self.assertEqual(
                    run_command("uma_alignment", [str(evidence)], worker), 0
                )
            values = json.loads(evidence.read_text())
            expected = Path(values["managed_temp"])
            self.assertEqual(Path(values["python_temp"]), expected)
            self.assertEqual(Path(values["java_temp"]), expected)
            self.assertEqual(Path(values["java_file"]).parent, expected)
            self.assertFalse(expected.exists())
            record = json.loads(
                next(runtime_home.glob("runs/*/run.json")).read_text()
            )
            self.assertEqual(record["cleanup"]["status"], "CLEANED")

    def test_five_commands_handoff_partial_results_inside_uma_assay(self):
        import numpy as np
        import openpyxl
        import tifffile

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "images"
            source.mkdir()
            runtime_home = Path(temporary) / "runtime"
            environment = dict(os.environ, UMA_RUNTIME_HOME=str(runtime_home))
            stack = np.zeros((7, 32, 32), dtype=np.uint16)
            stack[1:6, 5:27, 8:12] = 30000
            stack[2:5, 9:23, 21:25] = 45000
            tifffile.imwrite(
                source / "good_WellB02.tiff",
                stack,
                imagej=True,
                resolution=(2, 2),
                metadata={"axes": "ZYX", "spacing": 1.0, "unit": "um"},
            )
            (source / "bad_WellB03.tiff").write_bytes(b"unreadable TIFF")
            (source / "._good_WellB02.tiff").write_bytes(b"metadata")
            config = source / "input.json"
            config.write_text(json.dumps({"folder_paths": [str(source)]}))
            for command, answers, pattern in (
                (
                    "uma_alignment",
                    "1\ny\n",
                    "uma_assay/Alignment_assay_results_*",
                ),
                (
                    "uma_thickness",
                    "2\n1\ny\n",
                    "uma_assay/Thickness_assay_results_*",
                ),
                ("area_analysis", "1\n", "uma_assay/Area_assay_results_*"),
            ):
                with self.subTest(command=command):
                    result = subprocess.run(
                        [
                            str(Path(sys.executable).parent / command),
                            "-i",
                            str(config),
                        ],
                        input=answers,
                        text=True,
                        env=environment,
                        capture_output=True,
                        timeout=300,
                    )
                    self.assertEqual(
                        result.returncode, 1, result.stdout + result.stderr
                    )
                    output = next(source.glob(pattern))
                    status = json.loads(
                        (output / "run_status.json").read_text()
                    )
                    self.assertEqual(
                        status["status"],
                        "PARTIAL",
                        result.stdout + result.stderr,
                    )
                    self.assertEqual(status["processed_images"], 1)
                    self.assertEqual(status["failed_images"], 1)
                    self.assertEqual(status["unprocessed_images"], 0)
                    errors = (output / "image_errors.csv").read_text()
                    self.assertIn("bad_WellB03.tiff", errors)
                    self.assertNotIn("._good", errors)
                    self.assertTrue(status["summary_sha256"])
                    if command == "uma_alignment":
                        journal = (
                            source / "uma_assay/UMA_Logs/1_alignment.log"
                        ).read_text()
                        for number, stage, count in (
                            (1, "Projection", 2),
                            (2, "Orientation", 1),
                            (3, "Summary", 1),
                        ):
                            self.assertIn(
                                f"Stage {number}/3 | {stage}: "
                                f"{count}/{count} finished",
                                journal,
                            )
                    if command == "uma_thickness":
                        journal = (
                            source / "uma_assay/UMA_Logs/2_thickness.log"
                        ).read_text()
                        for number, operation in enumerate(
                            (
                                "Opening image...",
                                "Extracting fibronectin channel...",
                                "Performing Reslice...",
                                "Performing Z projection...",
                                "Applying Maximum filter...",
                                "Applying Gaussian Blur...",
                                "Subtracting background...",
                                "Applying threshold...",
                                "Running Local Thickness...",
                                "Measuring thickness...",
                                "Closed all images.",
                            ),
                            1,
                        ):
                            self.assertRegex(
                                journal,
                                r"Stage 1/1 \| Thickness: [01]/2 finished"
                                r"(?: \| 1 failed)? \| "
                                + f"Operation {number}/11: "
                                + re.escape(operation)
                                + r" \| good_WellB02\.tiff",
                            )
                        self.assertIn(
                            "Thickness: 2/2 finished | 1 failed", journal
                        )

            container = source / "uma_assay"
            collected = subprocess.run(
                [
                    str(Path(sys.executable).parent / "uma_collect_results"),
                    "-i",
                    str(config),
                ],
                text=True,
                env=environment,
                capture_output=True,
                timeout=60,
            )
            self.assertEqual(
                collected.returncode, 0, collected.stdout + collected.stderr
            )
            combined = next(container.glob("Combined_Results_*"))
            collection = json.loads((combined / "run_status.json").read_text())
            self.assertEqual(collection["retained_images"], 1)
            self.assertEqual(collection["excluded_images"], 1)
            self.assertEqual(collection["source_folder"], str(source))

            workbook = openpyxl.Workbook()
            sheet = workbook.active
            for column in range(1, 13):
                sheet.cell(1, column + 1, column)
            for row, label in enumerate("ABCDEFGH", 2):
                sheet.cell(row, 1, label)
            sheet["C3"] = "Control"
            workbook.save(combined / "Any plate filename.xlsx")
            workbook.close()

            reported = subprocess.run(
                [
                    str(Path(sys.executable).parent / "uma_report"),
                    "-i",
                    str(config),
                    "--fn-threshold",
                    "0",
                ],
                text=True,
                env=environment,
                capture_output=True,
                timeout=180,
            )
            self.assertEqual(
                reported.returncode, 0, reported.stdout + reported.stderr
            )
            output = next(container.glob("UMA_Report_*"))
            report = json.loads((output / "run_status.json").read_text())
            self.assertEqual(report["status"], "SUCCESS")
            self.assertEqual(report["source_folder"], str(source))
            self.assertEqual(report["combined_results_folder"], str(combined))
            self.assertEqual(report["total_images"], 1)
            self.assertEqual(report["processing_excluded_images"], 1)
            self.assertEqual(report["generated_plots"], 14)
            self.assertIn(f"UMA_Report_{source.name}_", output.name)
            self.assertEqual(list(combined.glob("UMA_Report_*")), [])
            self.assertEqual(
                {path.name for path in (container / "UMA_Logs").glob("*.log")},
                {
                    "1_alignment.log",
                    "2_thickness.log",
                    "3_area.log",
                    "4_collect_results.log",
                    "5_report.log",
                },
            )
            self.assertEqual(
                {path.name for path in source.iterdir() if path.is_dir()},
                {"uma_assay"},
            )
            records = [
                json.loads(path.read_text())
                for path in runtime_home.glob("runs/*/run.json")
            ]
            self.assertEqual(len(records), 5)
            self.assertEqual(
                sorted(record["exit_code"] for record in records),
                [0, 0, 1, 1, 1],
            )
            for record in records:
                with self.subTest(runtime_command=record["command"]):
                    self.assertTrue(record["completion"]["normal_completion"])
                    self.assertEqual(record["cleanup"]["status"], "CLEANED")
                    for removed in record["cleanup"]["removed"]:
                        self.assertFalse(Path(removed).exists())
            self.assertEqual(list((runtime_home / "tmp").iterdir()), [])
            self.assertEqual(list(source.rglob(".uma_tmp_*")), [])
