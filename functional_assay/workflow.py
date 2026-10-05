"""Output locations, journals and live progress for functional assays."""

from __future__ import annotations

import logging
import sys
import traceback
from contextvars import ContextVar
from pathlib import Path

from uma_tools.files import save_csv
from uma_tools.progress import CompactProgress
from uma_tools.run import RunLog as BaseRunLog
from uma_tools.run import utc_now

LOG_NAMES = {
    "stitching": "1_stitching.log",
    "cell_count": "2_cell_count.log",
    "functional_report": "3_functional_report.log",
    "survival_report": "4_survival_report.log",
}
EXCLUSION_COLUMNS = ("Well", "Stage", "Reason")
_ACTIVE_LOG = ContextVar("functional_run_log", default=None)


class NoInputError(ValueError):
    """The selected source contains no eligible inputs."""


def assay_directory(source: Path, *, create: bool = True) -> Path:
    """Use only the new layout; never create a missing image source."""
    if not source.is_dir():
        raise NoInputError(f"Source folder not found: {source}")
    directory = source / "uma_functional_assay"
    if directory.is_symlink():
        raise ValueError(
            f"Results directory must not be a symlink: {directory}"
        )
    if create:
        directory.mkdir(exist_ok=True)
    return directory


def _journal(source: Path, step: str):
    directory = assay_directory(source) / "UMA_Logs"
    if directory.is_symlink():
        raise ValueError(f"Log directory must not be a symlink: {directory}")
    directory.mkdir(exist_ok=True)
    path = directory / LOG_NAMES[step]
    if path.is_symlink():
        raise ValueError(f"Log file must not be a symlink: {path}")
    if path.exists():
        archive = directory / "archive"
        if archive.is_symlink():
            raise ValueError(f"Log archive must not be a symlink: {archive}")
        archive.mkdir(exist_ok=True)
        stamp = utc_now().replace(":", "").replace("+", "_")
        destination = archive / f"{path.stem}_{stamp}.log"
        # Never overwrite an existing archived journal.
        with destination.open("x", encoding="utf-8") as stream:
            stream.write(path.read_text(encoding="utf-8"))
    return path.open("w", encoding="utf-8", buffering=1)


class _PythonLogHandler(logging.Handler):
    def __init__(self, owner):
        super().__init__()
        self.owner = owner

    def emit(self, record):
        self.owner.event(
            record.levelname,
            record.name,
            self.format(record),
            console=record.levelno >= logging.WARNING,
        )


class RunLog(BaseRunLog):
    """Mirror run events to a rotating journal and one live terminal line."""

    def __init__(self, output: Path, source: Path, step: str):
        self.source = source
        self.journal = _journal(source, step)
        try:
            super().__init__(output)
        except BaseException:
            self.journal.close()
            raise
        self.progress = CompactProgress()
        self.progress.__enter__()
        self.phase_label = source.name
        self.closed = False
        self.token = _ACTIVE_LOG.set(self)
        self.handler = _PythonLogHandler(self)
        self.root = logging.getLogger()
        self.previous_handlers = self.root.handlers[:]
        self.previous_level = self.root.level
        self.root.handlers = [self.handler]
        self.root.setLevel(logging.INFO)

    def event(self, level, stage, message, console=True, timestamp=None):
        row = super().event(level, stage, message, False, timestamp)
        line = f"[{row['Timestamp_UTC']}] [{level}] [{stage}] {message}"
        self.journal.write(line + "\n")
        self.journal.flush()
        if console and (level != "INFO" or stage in {"Output", "Plan"}):
            self.progress.message(f"{level}: {stage}: {message}")
        return row

    def phase(
        self,
        index,
        total,
        name,
        *,
        finished=None,
        count=None,
        unit="wells",
        detail="",
    ):
        label = f"{self.source.name}: stage {index}/{total} {name}"
        if finished is not None:
            label += f" | {finished}/{count} {unit} finished"
        self.phase_label = label
        self.activity(detail)
        self.event("INFO", "Progress", self.progress.label, console=False)
        if not self.progress.tty:
            self.progress.message(self.progress.label)

    def activity(self, detail):
        self.progress.update(
            self.phase_label + (f" | {detail}" if detail else "")
        )

    def record_error(self, stage, error):
        self.event("ERROR", stage, str(error) or type(error).__name__)
        self.event("ERROR", "Traceback", traceback.format_exc(), console=False)

    def close(self):
        if self.closed:
            return
        self.closed = True
        self.root.handlers = self.previous_handlers
        self.root.setLevel(self.previous_level)
        _ACTIVE_LOG.reset(self.token)
        self.progress.__exit__(None, None, None)
        try:
            super().close()
        finally:
            self.journal.close()


def activity(detail: str) -> None:
    """Update the active well/operation without resetting its counter."""
    log = _ACTIVE_LOG.get()
    if log is not None:
        log.activity(detail)


def warning(message: str) -> None:
    log = _ACTIVE_LOG.get()
    if log is not None:
        log.event("WARNING", "Input discovery", message)
    else:
        print(f"WARNING: {message}", flush=True)


def plot_progress(finished: int, total: int, detail: str = "") -> None:
    log = _ACTIVE_LOG.get()
    if log is not None:
        log.phase(
            3,
            4,
            "Plots",
            finished=finished,
            count=total,
            unit="plots",
            detail=detail,
        )


def exclusions(status: dict) -> list[dict]:
    """Use recorded exclusions, never infer failed wells from zero values."""
    return [
        {
            "Well": well,
            "Stage": row.get("stage", "Cell analysis"),
            "Reason": row.get("error", "Unsuccessful well"),
        }
        for well, row in status.get("wells", {}).items()
        if row.get("status") != "completed"
    ]


def save_exclusions(output: Path, rows: list[dict]) -> None:
    save_csv(
        output / "Processing_Exclusions.csv",
        EXCLUSION_COLUMNS,
        rows,
        encoding="utf-8-sig",
    )


def outcome(completed: int, failures: int) -> str:
    return (
        "PARTIAL"
        if completed and failures
        else ("SUCCESS" if completed else "FAILED")
    )


def command_error(step: str, error: BaseException, source: Path | None = None):
    """Retain startup/output/shutdown failures even without a run folder."""
    state = (
        "CANCELLED"
        if isinstance(error, KeyboardInterrupt)
        else ("NO_INPUT" if isinstance(error, NoInputError) else "ERROR")
    )
    message = str(error) or type(error).__name__
    print(f"{state}: {step}: {message}", file=sys.stderr, flush=True)
    destination = (
        source if source is not None and source.is_dir() else Path.cwd()
    )
    try:
        directory = assay_directory(destination) / "UMA_Logs"
        directory.mkdir(exist_ok=True)
        path = directory / LOG_NAMES[step]
        if directory.is_symlink() or path.is_symlink():
            raise ValueError("Refusing to follow a log symlink")
        with path.open("a", encoding="utf-8") as stream:
            stream.write(f"[{utc_now()}] [{state}] {message}\n")
            stream.write(traceback.format_exc() + "\n")
        print(f"Diagnostics: {path}", file=sys.stderr, flush=True)
    except (OSError, ValueError) as log_error:
        print(f"Could not save diagnostics: {log_error}", file=sys.stderr)


def batch_error(step: str, error: BaseException, input_path: str):
    """Put command-boundary failures in each available source journal."""
    from uma_tools.config import read_config

    try:
        folders = read_config(Path(input_path))
    except (OSError, ValueError):
        folders = []
    available = list(
        dict.fromkeys(folder for folder in folders if folder.is_dir())
    )
    for folder in available or [None]:
        command_error(step, error, folder)
