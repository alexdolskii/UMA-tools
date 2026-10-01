"""Native command recovery and process exit, isolated from the test JVM."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


@unittest.skipUnless(
    os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1", "Requires the Fiji runtime"
)
class PartialFijiTests(unittest.TestCase):
    def test_three_commands_keep_good_images_after_an_unreadable_tiff(self):
        import numpy as np
        import tifffile

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
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
                ("uma_alignment", "1\ny\n", "Alignment_assay_results_*"),
                ("uma_thickness", "2\n1\ny\n", "Thickness_assay_results_*"),
                ("area_analysis", "1\n", "Area_assay_results_*"),
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
