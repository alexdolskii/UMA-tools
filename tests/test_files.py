"""Direct unit tests for the shared file serialization primitives."""

import csv as csv_module
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from uma_tools import files


class SaveJsonTests(unittest.TestCase):
    def test_atomic_default_writes_indented_utf8_with_a_trailing_newline(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.json"
            files.save_json(path, {"b": 1, "a": 2})
            text = path.read_text(encoding="utf-8")
            self.assertTrue(text.endswith("\n"))
            self.assertEqual(json.loads(text), {"b": 1, "a": 2})
            self.assertEqual(list(path.parent.glob("*.tmp")), [])

    def test_can_disable_the_atomic_rename_and_trailing_newline(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.json"
            files.save_json(
                path, {"a": 1}, atomic=False, trailing_newline=False
            )
            text = path.read_text(encoding="utf-8")
            self.assertFalse(text.endswith("\n"))

    def test_rejects_nonfinite_values_unless_allow_nan_is_set(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.json"
            with self.assertRaises(ValueError):
                files.save_json(path, {"a": float("nan")}, allow_nan=False)
            files.save_json(path, {"a": float("nan")})
            self.assertIn("NaN", path.read_text(encoding="utf-8"))

    def test_leaves_a_recognizable_temporary_file_if_the_rename_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.json"

            def fail_rename(self_path, target):
                raise OSError("deliberate rename failure")

            with patch.object(Path, "replace", fail_rename):
                with self.assertRaises(OSError):
                    files.save_json(
                        path, {"a": 1}, temporary_suffix=".pending"
                    )
            self.assertFalse(path.exists())
            self.assertTrue(path.with_name("data.json.pending").exists())


class SaveCsvTests(unittest.TestCase):
    def test_writes_the_header_and_rows_in_the_requested_column_order(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.csv"
            files.save_csv(
                path,
                ["B", "A"],
                [{"A": "1", "B": "2"}, {"A": "3", "B": "4"}],
            )
            with path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv_module.reader(stream))
            self.assertEqual(rows, [["B", "A"], ["2", "1"], ["4", "3"]])

    def test_supports_utf8_sig_encoding_for_excel_compatibility(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.csv"
            files.save_csv(
                path, ["Name"], [{"Name": "тест"}], encoding="utf-8-sig"
            )
            self.assertTrue(path.read_bytes().startswith(b"\xef\xbb\xbf"))


class Sha256FileTests(unittest.TestCase):
    def test_matches_hashlib_across_multiple_read_chunks(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.bin"
            content = os.urandom(1024 * 1024 + 17)
            path.write_bytes(content)
            self.assertEqual(
                files.sha256_file(path), hashlib.sha256(content).hexdigest()
            )


class SafeLabelTests(unittest.TestCase):
    def test_replaces_unsafe_characters_and_trims_dots_and_spaces(self):
        self.assertEqual(
            files.safe_label('a/b\\c:d*e?f"g<h>i|j'), "a_b_c_d_e_f_g_h_i_j"
        )
        self.assertEqual(files.safe_label("  name.  "), "name")

    def test_falls_back_to_images_for_empty_or_dot_only_names(self):
        self.assertEqual(files.safe_label(""), "images")
        self.assertEqual(files.safe_label("..."), "images")

    def test_shortens_long_names_with_a_deterministic_hash_suffix(self):
        long_name = "x" * 200
        label = files.safe_label(long_name)
        self.assertLessEqual(len(label.encode("utf-8")), 109)
        self.assertEqual(label, files.safe_label(long_name))


if __name__ == "__main__":
    unittest.main()
