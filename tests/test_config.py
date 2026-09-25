"""Direct unit tests for the shared configuration primitives."""

import json
import os
import tempfile
import unittest
from pathlib import Path

from uma_tools import config


class RejectMetadataJsonTests(unittest.TestCase):
    def test_rejects_appledouble_names_without_reading_the_file(self):
        # The path does not exist; success without a FileNotFoundError
        # proves the check is a pure filename check.
        path = Path("/nonexistent/._input.json")
        with self.assertRaisesRegex(config.ValidationError, "metadata"):
            config.reject_metadata_json(path)

    def test_allows_regular_names(self):
        config.reject_metadata_json(Path("/nonexistent/input.json"))

    def test_supports_a_custom_error_type_and_message(self):
        class CustomError(Exception):
            pass

        with self.assertRaisesRegex(CustomError, "custom message"):
            config.reject_metadata_json(
                Path("._x.json"),
                error_type=CustomError,
                message="custom message",
            )


class LoadJsonTests(unittest.TestCase):
    def test_reads_utf8_sig_bom_and_parses_data(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "input.json"
            path.write_bytes(
                b"\xef\xbb\xbf" + json.dumps({"a": 1}).encode("utf-8")
            )
            self.assertEqual(config.load_json(path), {"a": 1})

    def test_rejects_metadata_before_reading_by_default(self):
        path = Path("/nonexistent/._input.json")
        with self.assertRaisesRegex(config.ValidationError, "metadata"):
            config.load_json(path)

    def test_can_disable_metadata_rejection(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "._input.json"
            path.write_text(json.dumps({"a": 1}), encoding="utf-8")
            self.assertEqual(
                config.load_json(path, reject_metadata=False), {"a": 1}
            )


class FolderEntriesTests(unittest.TestCase):
    def test_requires_a_nonempty_list(self):
        for data in ({"folder_paths": []}, {}, None, {"folder_paths": "x"}):
            with self.subTest(data=data):
                with self.assertRaises(config.ValidationError):
                    config.folder_entries(data)

    def test_rejects_blank_or_non_string_entries(self):
        for bad in ([""], ["   "], [123], [None]):
            with self.subTest(bad=bad):
                with self.assertRaises(config.ValidationError):
                    config.folder_entries({"folder_paths": bad})

    def test_preserves_order_and_duplicates(self):
        values = ["b", "a", "a"]
        self.assertEqual(
            config.folder_entries({"folder_paths": values}), values
        )


class ResolvePathTests(unittest.TestCase):
    def test_keeps_absolute_paths_only_normalized(self):
        self.assertEqual(
            config.resolve_path("/tmp/some/../folder", Path("/unused")),
            Path("/tmp/some/../folder").resolve(),
        )

    def test_joins_relative_paths_to_the_given_root(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "images").mkdir()
            self.assertEqual(
                config.resolve_path("images", root),
                (root / "images").resolve(),
            )

    def test_expands_the_user_home_directory(self):
        self.assertEqual(
            config.resolve_path("~", Path("/unused")), Path.home().resolve()
        )


class ResolveFoldersTests(unittest.TestCase):
    def test_preserves_order_and_repeats_without_following_symlinks(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            target = root / "images"
            target.mkdir()
            link = root / "link"
            link.symlink_to(target, target_is_directory=True)
            resolved = config.resolve_folders(
                ["link", "link"], relative_to=root
            )
            self.assertEqual(resolved, [root / "link", root / "link"])

    def test_can_follow_symlinks_when_requested(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            target = root / "images"
            target.mkdir()
            link = root / "link"
            link.symlink_to(target, target_is_directory=True)
            resolved = config.resolve_folders(
                ["link"], relative_to=root, use_realpath=True
            )
            self.assertEqual(resolved, [target.resolve()])

    def test_defaults_to_the_current_working_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            previous = Path.cwd()
            try:
                os.chdir(root)
                self.assertEqual(
                    config.resolve_folders(["images"]),
                    [Path(os.path.abspath("images"))],
                )
            finally:
                os.chdir(previous)


class ReadConfigTests(unittest.TestCase):
    def test_resolves_folders_relative_to_the_json_argument(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            images = root / "images"
            images.mkdir()
            manifest = root / "input.json"
            manifest.write_text(json.dumps({"folder_paths": ["images"]}))
            previous = Path.cwd()
            try:
                os.chdir(root)
                self.assertEqual(
                    config.read_config(manifest), [images.resolve()]
                )
            finally:
                os.chdir(previous)

    def test_rejects_a_non_json_suffix(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "input.txt"
            path.write_text("{}")
            with self.assertRaises(config.ValidationError):
                config.read_config(path)

    def test_rejects_metadata_before_reading(self):
        path = Path("/nonexistent/._input.json")
        with self.assertRaisesRegex(config.ValidationError, "metadata"):
            config.read_config(path)


if __name__ == "__main__":
    unittest.main()
