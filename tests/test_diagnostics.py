"""Diagnostics inspect bounded known locations without removing data."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from uma_tools import diagnostics, runtime


class DiagnosticsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.home = self.root / "runtime"
        self.source = self.root / "images"
        self.assay = self.source / "uma_assay"
        self.assay.mkdir(parents=True)
        self.cache = self.root / "cache"
        self.cache.mkdir()
        patch = mock.patch.dict(os.environ, UMA_RUNTIME_HOME=str(self.home))
        patch.start()
        self.addCleanup(patch.stop)

    def test_diagnostics_and_runtime_imports_are_lightweight(self):
        script = (
            "import sys; from uma_tools import diagnostics, runtime; "
            "assert not {'imagej','jpype','scyjava','numpy','pandas',"
            "'openpyxl','matplotlib','numba'} & set(sys.modules)"
        )
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=15,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_read_only_scan_copies_and_archives_logs_using_folder_paths(self):
        protected = {
            self.cache / "component.jar": b"permanent dependency",
            self.assay / "summary.partial.csv": b"scientific checkpoint",
            self.assay / "foreign.tmp": b"unattributed data",
        }
        for path, content in protected.items():
            path.write_bytes(content)
        config = self.root / "input.json"
        config.write_text(json.dumps({"folder_paths": [str(self.source)]}))
        with (
            mock.patch.object(
                diagnostics, "cache_locations", return_value=[self.cache]
            ),
            mock.patch.object(
                diagnostics,
                "system_temporaries",
                return_value={"files": [], "checked_locations": []},
            ),
        ):
            self.assertIn(diagnostics.main(["-i", str(config)]), (0, 1))
            first = (self.home / "uma_diagnostics.log").read_bytes()
            copied = self.assay / "UMA_Logs" / "uma_diagnostics.log"
            self.assertEqual(first, copied.read_bytes())
            self.assertIn(b"foreign.tmp", first)
            self.assertNotIn(b"summary.partial.csv", first)
            self.assertIn(b"PERSISTENT_CACHE", first)
            self.assertIn(diagnostics.main(["-i", str(config)]), (0, 1))
        for folder in (self.home, self.assay / "UMA_Logs"):
            archives = list((folder / "archive").glob("uma_diagnostics_*.log"))
            self.assertEqual(len(archives), 1)
            self.assertEqual(archives[0].read_bytes(), first)
        for path, content in protected.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertFalse((self.home / "runs").exists())

    def test_bounded_scan_does_not_follow_symlinks(self):
        for index in range(5):
            (self.cache / f"{index}.bin").write_bytes(b"123")
        (self.assay / "linked_cache").symlink_to(
            self.cache, target_is_directory=True
        )
        errors = []
        self.assertEqual(
            diagnostics.size_of(self.assay, 20, errors),
            {"files": 0, "bytes": 0},
        )
        limited = diagnostics.size_of(self.cache, 2, errors)
        self.assertEqual(limited, {"files": 2, "bytes": 6})
        self.assertIn("Scan limit", errors[-1])
        self.assertEqual(len(list(self.cache.iterdir())), 5)

    def test_corrupt_record_is_reported_as_unknown_without_hiding_other_runs(
        self,
    ):
        run_id = "20261001T010203123456Z_abcdef123456"
        directory = self.home / "runs" / run_id
        directory.mkdir(parents=True)
        (directory / "run.json").write_text(
            json.dumps(
                {
                    "schema": runtime.SCHEMA,
                    "run_id": run_id,
                    "processes": ["invalid process"],
                }
            )
        )
        errors = []
        runs = diagnostics.managed_runs(self.home, 100, errors)
        self.assertEqual(runs[0]["status"], "UNKNOWN")
        self.assertIn("Invalid process registry", errors[0])
        self.assertTrue((directory / "run.json").is_file())

    def test_metadata_json_and_relative_runtime_home_are_rejected(self):
        with self.assertRaises(SystemExit) as error:
            diagnostics.main(["-i", str(self.root / "._input.json")])
        self.assertEqual(error.exception.code, 2)
        with mock.patch.dict(os.environ, UMA_RUNTIME_HOME="relative"):
            with self.assertRaisesRegex(ValueError, "absolute path"):
                runtime.runtime_home()
        self.assertFalse(self.home.exists())

    def test_linked_log_destination_is_not_overwritten(self):
        self.home.mkdir()
        protected = self.root / "important.txt"
        protected.write_text("keep")
        (self.home / "uma_diagnostics.log").symlink_to(protected)
        with self.assertRaisesRegex(OSError, "linked diagnostic log"):
            diagnostics.save_log(self.home, {})
        self.assertEqual(protected.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
