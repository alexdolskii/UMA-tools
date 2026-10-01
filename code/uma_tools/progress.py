"""Command journals, compact progress and recoverable terminal prompts.

The context is active only for the five core commands. Shared library
callers retain their existing logging and console behavior.
"""

from __future__ import annotations

import functools
import inspect
import logging
import shutil
import sys
import threading
import time
import traceback
from contextlib import contextmanager
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path

from . import package_version
from .files import assay_directory
from .runtime import runtime_context

_SESSION = ContextVar("uma_session", default=None)
_FOLDER = ContextVar("uma_folder", default=None)
_IMAGE_PROGRESS = ContextVar("uma_image_progress", default=None)
STEPS = {
    "alignment": "1_alignment.log",
    "thickness": "2_thickness.log",
    "area": "3_area.log",
    "collect_results": "4_collect_results.log",
    "report": "5_report.log",
}


class Cancelled(Exception):
    """
    A deliberate terminal cancellation, rather than an analysis error.
    """


CANCELLATIONS = (Cancelled, KeyboardInterrupt, EOFError)


class CompactProgress:
    """
    Render elapsed activity without inventing a completion percentage.
    """

    def __init__(self, stream=None):
        self.stream = stream if stream is not None else sys.stdout
        self.tty = self.stream.isatty()
        self.label = "Starting"
        self.started = time.monotonic()
        self.lock = threading.RLock()
        self.stop = threading.Event()
        self.thread = None
        self.paused = False

    def __enter__(self):
        if self.tty:
            self.thread = threading.Thread(target=self._animate, daemon=True)
            self.thread.start()
        return self

    def _animate(self):
        while not self.stop.wait(0.2):
            self.render()

    def render(self):
        with self.lock:
            if self.tty and not self.paused:
                elapsed = int(time.monotonic() - self.started)
                text = f"{self.label} | elapsed {elapsed}s"
                width = max(20, shutil.get_terminal_size((80, 24)).columns - 1)
                self.stream.write("\r\033[2K" + text[:width])
                self.stream.flush()

    def update(self, label):
        with self.lock:
            self.label = str(label).replace("\n", " ")
            self.render()

    def message(self, message):
        with self.lock:
            if self.tty:
                self.stream.write("\r\033[2K")
            self.stream.write(str(message) + "\n")
            self.stream.flush()
            self.render()

    def __exit__(self, *_):
        self.stop.set()
        if self.thread is not None:
            self.thread.join(timeout=1)
        if self.tty:
            with self.lock:
                self.stream.write("\r\033[2K")
                self.stream.flush()


class JournalHandler(logging.Handler):
    """
    Forward Python records, including tracebacks, to the active journal.
    """

    def emit(self, record):
        session = current_session()
        if session is not None:
            session.event(
                record.levelname,
                record.name,
                self.format(record),
                announce=False,
            )
            if record.levelno >= logging.WARNING:
                session.progress.message(
                    f"{record.levelname}: {record.getMessage()}"
                )


class CommandSession:
    """
    Own exactly one current log per source folder for this command.
    """

    def __init__(self, step, folders=(), input_json=""):
        self.step = step
        self.folders = list(dict.fromkeys(Path(p).resolve() for p in folders))
        self.input_json = str(input_json)
        self.streams = {}
        self.unavailable = {}
        self.outcomes = {}
        self.progress = CompactProgress()
        self.exit_code = 1
        self.token = None
        self.root_handlers = []
        self.root_level = logging.WARNING

    def __enter__(self):
        # Start animation on the first processing phase, after prompts.

        self.token = _SESSION.set(self)
        root = logging.getLogger()
        self.root_handlers = root.handlers[:]
        self.root_level = root.level
        root.handlers = [JournalHandler()]
        root.setLevel(logging.INFO)
        try:
            if self.folders:
                self.register(self.folders, self.input_json)
        except BaseException:
            self.__exit__(*sys.exc_info())
            raise
        return self

    def register(self, folders, input_json):
        self.folders = list(dict.fromkeys(Path(p).resolve() for p in folders))
        self.input_json = str(input_json)
        for folder in self.folders:
            try:
                if not folder.is_dir():
                    raise FileNotFoundError(
                        f"Source folder not found: {folder}"
                    )
                logs = assay_directory(folder, create=True) / "UMA_Logs"
                logs.mkdir(exist_ok=True)
                path = logs / STEPS[self.step]
                if path.exists():
                    archive = logs / "archive"
                    archive.mkdir(exist_ok=True)
                    stamp = datetime.now(timezone.utc).strftime(
                        "%Y%m%d_%H%M%S_%f"
                    )
                    target = archive / f"{path.stem}_{stamp}.log"
                    counter = 0
                    while target.exists():
                        counter += 1
                        target = archive / f"{path.stem}_{stamp}_{counter}.log"
                    path.rename(target)
                self.streams[folder] = path.open(
                    "x", encoding="utf-8", buffering=1
                )
            except OSError as error:
                self.unavailable[folder] = str(error)
                self.progress.message(f"ERROR: {error}")
        self.event(
            "STARTED",
            self.step,
            f"UMA-tools {package_version()}; input={self.input_json}; "
            f"Python={sys.version.split()[0]}; executable={sys.executable}",
        )
        context = runtime_context()
        if context:
            self.event("INFO", "Temporary resources", context)
        self.progress.message(
            f"UMA {self.step}: {len(self.folders)} source folder(s). "
            "Detailed logs: <source>/uma_assay/UMA_Logs/" + STEPS[self.step]
        )

    def event(self, level, stage, message, *, announce=False, folder=None):
        target = Path(folder).resolve() if folder else _FOLDER.get()
        streams = (
            [self.streams[target]]
            if target in self.streams
            else list(self.streams.values())
            if target is None
            else []
        )
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        line = f"[{stamp}] [{level}] [{stage}] {message}"
        for stream in streams:
            try:
                stream.write(line + "\n")
                stream.flush()
            except OSError as error:
                for source, owned in list(self.streams.items()):
                    if owned is stream:
                        self.unavailable[source] = (
                            f"Journal write failed: {error}"
                        )
                        del self.streams[source]
                try:
                    stream.close()
                except OSError:
                    pass
                raise
        if announce:
            self.progress.message(f"{level}: {message}")

    def __exit__(self, exc_type, exc, tb):
        try:
            if exc_type is not None:
                cancelled = issubclass(exc_type, CANCELLATIONS)
                self.exit_code = 130 if cancelled else 1
                self.event(
                    "CANCELLED" if cancelled else "FAILED",
                    self.step,
                    "Cancelled by user"
                    if cancelled
                    else "".join(
                        traceback.format_exception(exc_type, exc, tb)
                    ),
                )
            self.event(
                "FINISHED",
                self.step,
                f"Command exit code: {self.exit_code}",
            )
        finally:
            for stream in self.streams.values():
                stream.close()
            root = logging.getLogger()
            for handler in root.handlers:
                if handler not in self.root_handlers:
                    handler.close()
            root.handlers = self.root_handlers
            root.setLevel(self.root_level)
            _SESSION.reset(self.token)
            self.progress.__exit__()


def current_session():
    return _SESSION.get()


def current_folder():
    return _FOLDER.get()


@contextmanager
def folder_scope(folder):
    """
    Route messages to one source, and restore routing on every exit.
    """
    folder = Path(folder).resolve()
    previous = _FOLDER.get()
    token = _FOLDER.set(folder)
    session = current_session()
    try:
        if session is not None and previous != folder:
            if folder in session.unavailable:
                raise OSError(session.unavailable[folder])
            session.event("STARTED", "Folder", str(folder), announce=True)
            session.progress.update(folder.name)
        yield
    finally:
        _FOLDER.reset(token)


@contextmanager
def image_progress(counter, filename, *, operations=()):
    """Keep an image's counter visible while its operations change."""
    token = _IMAGE_PROGRESS.set(
        (current_folder(), counter, filename, tuple(operations))
    )
    try:
        yield
    finally:
        _IMAGE_PROGRESS.reset(token)


def update_activity(message, *, record=False):
    """Show an operation with any active image counter and filename."""
    session = current_session()
    folder = current_folder()
    image = _IMAGE_PROGRESS.get()
    if image is not None and image[0] == folder:
        message = str(message).strip()
        operations = image[3]
        if message in operations:
            index = operations.index(message) + 1
            message = f"Operation {index}/{len(operations)}: {message}"
        message = f"{image[1]} | {message} | {image[2]}"
    if session is None:
        if record:
            print(message, flush=True)
    else:
        if record:
            session.event("INFO", "Stage", message)
        label = f"{folder.name}: " if folder is not None else ""
        if session.progress.thread is None and session.progress.tty:
            session.progress.__enter__()
        session.progress.update(label + str(message).strip())


def phase(message):
    """Record an operation without flooding a redirected terminal."""
    update_activity(message, record=True)


def outcome(status, message):
    session = current_session()
    if session is None:
        print(f"{status}: {message}", flush=True)
    else:
        session.event(status, "Summary", message, announce=True)


def ask_choice(prompt, choices):
    """Retry a constrained answer; q, EOF and Ctrl-C cancel cleanly."""
    while True:
        value = prompt_value(prompt).strip().lower()
        if value in ("q", "quit"):
            raise Cancelled()
        if value in choices:
            return choices[value]
        print("Please enter " + ", ".join(choices) + " (q to cancel).")


def ask_channel():
    while True:
        value = prompt_value(
            "Enter fibronectin channel index (starting from 1): "
        ).strip()
        if value.lower() in ("q", "quit"):
            raise Cancelled()
        try:
            channel = int(value)
        except ValueError:
            channel = 0
        if channel >= 1:
            return channel
        print("Please enter a positive integer (q to cancel).")


def confirm_start():
    if not ask_choice(
        "Do you want to start processing? (y/n): ",
        {"y": True, "yes": True, "n": False, "no": False},
    ):
        raise Cancelled()


def register_sources(folders, input_json):
    session = current_session()
    if session is not None:
        session.register(folders, input_json)


def folder_logged(parameter):
    """
    Route a folder-level operation without changing its public API.
    """

    def decorate(function):
        signature = inspect.signature(function)

        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            bound = signature.bind(*args, **kwargs)
            folder = Path(bound.arguments[parameter]).resolve()
            session = current_session()
            if session is not None:
                session.outcomes[str(folder)] = "RUNNING"
            try:
                with folder_scope(folder):
                    result = function(*args, **kwargs)
                if session is not None:
                    session.outcomes[str(folder)] = getattr(
                        result, "value", result
                    )
                return result
            except BaseException as error:
                if session is not None:
                    session.outcomes[str(folder)] = (
                        "CANCELLED"
                        if isinstance(error, CANCELLATIONS)
                        else "ERROR"
                    )
                raise

        return wrapped

    return decorate


def prompt_value(prompt):
    session = current_session()
    if session is None:
        return input(prompt)
    with session.progress.lock:
        session.progress.paused = True
        if session.progress.tty:
            session.progress.stream.write("\r\033[2K")
    try:
        return input(prompt)
    finally:
        session.progress.paused = False


def folder_error(folder, error):
    session = current_session()
    if session is None:
        logging.getLogger(__name__).exception("Folder failed: %s", folder)
    else:
        session.event(
            "FAILED", "Folder", traceback.format_exc(), folder=folder
        )
        session.progress.message(f"FAILED: {folder}: {error}")


def console(*args, **kwargs):
    """
    Print a complete message without colliding with a TTY progress line.
    """
    session = current_session()
    if session is None:
        print(*args, **kwargs)
    else:
        text = kwargs.get("sep", " ").join(str(arg) for arg in args)
        session.progress.message(text)
