"""Direct unit tests for run directory allocation and diagnostics logging."""

import csv
import logging
import tempfile
import unittest
from pathlib import Path

from uma_tools import run


class UtcNowTests(unittest.TestCase):
    def test_returns_a_utc_offset_iso_timestamp(self):
        self.assertTrue(run.utc_now().endswith("+00:00"))

    def test_honors_the_requested_timespec(self):
        # "seconds" precision never includes a fractional-second component.
        value = run.utc_now(timespec="seconds")
        self.assertNotIn(".", value.split("+")[0])


class UniqueOutputTests(unittest.TestCase):
    def test_creates_a_new_directory_named_with_the_prefix_and_timestamp(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            run_id, output = run.unique_output(
                root, "Prefix_", timestamp="20260101_000000"
            )
            self.assertEqual(run_id, "20260101_000000")
            self.assertEqual(output, root / "Prefix_20260101_000000")
            self.assertTrue(output.is_dir())

    def test_appends_a_zero_padded_counter_when_the_run_id_collides(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "Prefix_20260101_000000").mkdir()
            run_id, output = run.unique_output(
                root, "Prefix_", timestamp="20260101_000000", counter_width=3
            )
            self.assertEqual(run_id, "20260101_000000_001")
            self.assertEqual(output, root / "Prefix_20260101_000000_001")

    def test_includes_the_process_id_when_requested(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            run_id, output = run.unique_output(
                root, "Prefix_", timestamp="20260101_000000", include_pid=True
            )
            self.assertTrue(run_id.startswith("20260101_000000_"))
            self.assertTrue(output.is_dir())

    def test_raises_the_requested_error_type_and_message_after_max_attempts(
        self,
    ):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "Prefix_20260101_000000").mkdir()
            with self.assertRaisesRegex(RuntimeError, "custom message"):
                run.unique_output(
                    root,
                    "Prefix_",
                    timestamp="20260101_000000",
                    max_attempts=1,
                    error_type=RuntimeError,
                    error_message="custom message",
                )

    def test_never_reuses_an_existing_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            first_id, first_output = run.unique_output(root, "Prefix_")
            second_id, second_output = run.unique_output(root, "Prefix_")
            self.assertNotEqual(first_output, second_output)
            self.assertNotEqual(first_id, second_id)


class LoggerTests(unittest.TestCase):
    def test_make_logger_writes_to_both_its_file_and_stdout_handlers(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            logger = run.make_logger(output)
            try:
                self.assertEqual(len(logger.handlers), 2)
                logger.info("hello")
                for handler in logger.handlers:
                    handler.flush()
                self.assertIn(
                    "hello", (output / "run.log").read_text(encoding="utf-8")
                )
            finally:
                run.close_logger(logger)
            self.assertEqual(logger.handlers, [])

    def test_scoped_file_log_removes_its_handler_even_after_an_error(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            logger = logging.getLogger(self.id())
            previous_level = logger.level
            try:
                with self.assertRaises(RuntimeError):
                    with run.scoped_file_log(logger, directory, "scoped.log"):
                        logger.info("inside the scope")
                        raise RuntimeError("boom")
                self.assertEqual(logger.handlers, [])
                self.assertEqual(logger.level, previous_level)
                self.assertIn(
                    "inside the scope",
                    (directory / "scoped.log").read_text(encoding="utf-8"),
                )
            finally:
                logger.setLevel(previous_level)


class RunLogTests(unittest.TestCase):
    def test_event_writes_matching_text_and_csv_records(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            log = run.RunLog(directory)
            try:
                row = log.event("INFO", "Stage", "message")
                self.assertEqual(row["Level"], "INFO")
                self.assertEqual(row["Stage"], "Stage")
                self.assertEqual(row["Message"], "message")
                self.assertEqual(log.events, [row])
            finally:
                log.close()
            text = (directory / "run.log").read_text(encoding="utf-8")
            self.assertIn("[INFO] [Stage] message", text)
            with (directory / "run_log.csv").open(
                encoding="utf-8-sig", newline=""
            ) as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["Message"], "message")

    def test_keep_events_false_still_writes_files_without_buffering_them(
        self,
    ):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            log = run.RunLog(directory, keep_events=False)
            try:
                log.event("INFO", "Stage", "message")
                self.assertEqual(log.events, [])
            finally:
                log.close()
            self.assertTrue((directory / "run.log").exists())

    def test_append_mode_does_not_repeat_the_csv_header(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            log = run.RunLog(directory)
            log.event("INFO", "Stage", "first")
            log.close()
            log = run.RunLog(directory, append=True)
            log.event("INFO", "Stage", "second")
            log.close()
            lines = (directory / "run_log.csv").read_text(
                encoding="utf-8-sig"
            ).splitlines()
            self.assertEqual(
                sum(1 for line in lines if line.startswith("Timestamp_UTC")),
                1,
            )
            self.assertEqual(len(lines), 3)


if __name__ == "__main__":
    unittest.main()
