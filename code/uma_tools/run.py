"""Run directories and scoped, flushed diagnostics for UMA workflows."""

from __future__ import annotations

import csv
import logging
import os
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .contracts import EVENT_COLUMNS


def utc_now(*, timespec: str = "auto") -> str:
    """Return an ISO UTC timestamp with the requested precision."""
    return datetime.now(timezone.utc).isoformat(timespec=timespec)


def unique_output(
    parent: Path,
    prefix: str,
    *,
    timestamp: str | None = None,
    include_pid: bool = False,
    counter_width: int = 0,
    max_attempts: int = 1000,
    error_type: type[Exception] = OSError,
    error_message: str = "Could not create a unique results directory",
) -> tuple[str, Path]:
    """Allocate a new run without creating a missing input folder.

    ``prefix`` includes its trailing underscore. Explicit timestamps
    retain the local-time naming policy of alignment and thickness;
    the other workflows use UTC. Existing folders are never reused.
    """
    stamp = timestamp or datetime.now(timezone.utc).strftime(
        "%Y%m%d_%H%M%S_%f"
    )
    base = f"{stamp}_{os.getpid()}" if include_pid else stamp
    for counter in range(max_attempts):
        suffix = f"_{counter:0{counter_width}d}" if counter else ""
        run_id = base + suffix
        output = parent / f"{prefix}{run_id}"
        try:
            output.mkdir()
            return run_id, output
        except FileExistsError:
            continue
    raise error_type(error_message)


def make_logger(output: Path) -> logging.Logger:
    """Create an isolated collector logger in its existing format."""
    logger = logging.Logger(str(output), level=logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (
        logging.FileHandler(output / "run.log", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def close_logger(logger: logging.Logger) -> None:
    """Close and remove all handlers owned by an isolated run logger."""
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)


@contextmanager
def scoped_file_log(
    logger: logging.Logger,
    directory: Path,
    filename: str = "log.log",
) -> Iterator[logging.FileHandler]:
    """Attach one folder's assay log and close it on every exit path."""
    handler = logging.FileHandler(directory / filename, mode="w")
    handler.setLevel(logging.INFO)
    handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    )
    previous_level = logger.level
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        yield handler
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous_level)


class RunLog:
    """Flush structured events to text, CSV and optionally the terminal.

    Area and report share this schema. Appending is used when recording
    ImageJ shutdown; reporting also retains the events for its workbook.
    """

    def __init__(
        self,
        directory: Path,
        append: bool = False,
        *,
        timespec: str = "seconds",
        keep_events: bool = True,
    ) -> None:
        self.events: list[dict[str, str]] = []
        self.keep_events = keep_events
        self.timespec = timespec
        self.path = directory / "run.log"
        mode = "a" if append else "w"
        csv_path = directory / "run_log.csv"
        needs_header = (
            not append or not csv_path.exists() or csv_path.stat().st_size == 0
        )
        self.text_stream = self.path.open(mode, encoding="utf-8", buffering=1)
        try:
            self.csv_stream = csv_path.open(
                mode, encoding="utf-8-sig", newline=""
            )
        except BaseException:
            self.text_stream.close()
            raise
        self.csv_writer = csv.DictWriter(
            self.csv_stream, fieldnames=EVENT_COLUMNS
        )
        self.writer = self.csv_writer
        if needs_header:
            self.csv_writer.writeheader()
        self.csv_stream.flush()

    def event(
        self,
        level: str,
        stage: str,
        message: object,
        console: bool = True,
        timestamp: str | None = None,
    ) -> dict[str, str]:
        """Persist an event and return its workbook record."""
        row = dict(
            zip(
                EVENT_COLUMNS,
                [
                    timestamp or utc_now(timespec=self.timespec),
                    level,
                    stage,
                    str(message),
                ],
            )
        )
        if self.keep_events:
            self.events.append(row)
        line = f"[{row['Timestamp_UTC']}] [{level}] [{stage}] {message}"
        self.text_stream.write(line + "\n")
        self.text_stream.flush()
        self.csv_writer.writerow(row)
        self.csv_stream.flush()
        if console:
            print(line, flush=True)
        return row

    def close(self) -> None:
        """Close both owned streams; repeated calls are harmless."""
        try:
            self.text_stream.close()
        finally:
            self.csv_stream.close()
