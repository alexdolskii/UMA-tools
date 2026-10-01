"""
Per-image evidence for recoverable failures and downstream exclusions.
"""

from __future__ import annotations

import logging
import traceback
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path

from . import package_version
from .files import save_csv, save_json, sha256_file
from .imagej import bioformats_log
from .progress import CANCELLATIONS, image_progress, outcome, phase
from .run import utc_now

SCHEMA = "uma-image-run-v1"
ERROR_COLUMNS = ["File_Name", "Stage", "Error"]
_ACTIVE = ContextVar("uma_image_run", default=None)
_LOG = logging.getLogger(__name__)


def image_names(folder, extensions):
    """
    Select only visible regular source images, never result directories.
    """
    return [
        item.name
        for item in Path(folder).iterdir()
        if not item.name.startswith(".")
        and item.is_file()
        and item.suffix.lower() in extensions
    ]


class ImageRun:
    """Publish a completion marker only after the summary is durable."""

    def __init__(self, output, source, names, parameters, status=None):
        self.output = Path(output)
        self.records = {
            name: {
                "File_Name": name,
                "Status": "PENDING",
                "Stage": "",
                "Error": "",
            }
            for name in names
        }
        self.status = status if status is not None else {}
        self.status.update(
            image_status_schema=SCHEMA,
            status="RUNNING",
            source_folder=str(source),
            package_version=package_version(),
            parameters=parameters,
            started_utc=utc_now(),
        )
        self.token = None
        self._stage_totals = {}
        self._stage_completed = {}
        _LOG.info("Analysis parameters: %s", parameters)
        self.save()

    def save(self):
        rows = list(self.records.values())
        self.status.update(
            images=rows,
            input_images=len(rows),
            processed_images=sum(r["Status"] == "SUCCESS" for r in rows),
            failed_images=sum(r["Status"] == "FAILED" for r in rows),
            unprocessed_images=sum(
                r["Status"] in ("PENDING", "RUNNING") for r in rows
            ),
        )
        save_json(self.output / "run_status.json", self.status)
        save_csv(
            self.output / "image_errors.csv",
            ERROR_COLUMNS,
            [
                {key: row[key] for key in ERROR_COLUMNS}
                for row in rows
                if row["Status"] == "FAILED"
            ],
        )

    def __enter__(self):
        self.token = _ACTIVE.set(self)
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            if exc_type is not None:
                self.status.update(
                    status="CANCELLED"
                    if isinstance(exc, CANCELLATIONS)
                    else "ERROR",
                    error=str(exc) or "Interrupted",
                    ended_utc=utc_now(),
                )
                self.save()
        finally:
            _ACTIVE.reset(self.token)

    def eligible(self, name):
        return name in self.records and (
            self.records[name]["Status"] != "FAILED"
        )

    def _progress_label(self, stage):
        """Count completed attempts in this stage, not final assay rows."""
        completed = self._stage_completed[stage]
        failed = sum(
            self.records[item]["Status"] == "FAILED" for item in completed
        )
        label = (
            f"{stage}: {len(completed)}/{self._stage_totals[stage]} finished"
        )
        if failed:
            label += f" | {failed} failed"
        return label

    def _show_progress(self, name, stage):
        phase(f"{self._progress_label(stage)} | {name}")

    @contextmanager
    def attempt(self, name, stage, *, final=False):
        if stage not in self._stage_totals:
            self._stage_totals[stage] = sum(
                self.eligible(item) for item in self.records
            )
            self._stage_completed[stage] = set()
        record = self.records[name]
        record.update(Status="RUNNING", Stage=stage)
        self.save()
        self._show_progress(name, stage)
        try:
            with (
                image_progress(self._progress_label(stage), name),
                bioformats_log(),
            ):
                yield
        except CANCELLATIONS:
            raise
        except Exception as error:
            record.update(Status="FAILED", Error=str(error))
            message = f"Image failed ({stage}): {name}: {error}"
            _LOG.exception(message)
            with (self.output / "image_tracebacks.log").open(
                "a", encoding="utf-8"
            ) as stream:
                stream.write(message + "\n" + traceback.format_exc() + "\n")
        else:
            if final:
                record["Status"] = "SUCCESS"
        self.save()
        self._stage_completed[stage].add(name)
        self._show_progress(name, stage)

    def finish(self, summary=None):
        self.save()
        if self.status["unprocessed_images"]:
            raise RuntimeError(
                "Some images have no terminal processing outcome"
            )
        count = self.status["processed_images"]
        failed = self.status["failed_images"]
        state = (
            "NO_INPUT"
            if not self.records
            else "PARTIAL"
            if count and failed
            else "SUCCESS"
            if count
            else "FAILED"
        )
        if summary is not None and count:
            self.status["summary_sha256"] = sha256_file(Path(summary))
        self.status.update(status=state, ended_utc=utc_now())
        self.save()
        outcome(
            state,
            f"{self.status['source_folder']}: {count} image(s) saved; "
            f"{failed} failed; {len(self.records)} selected. "
            f"Results: {self.output}",
        )
        return state

    def reject_collisions(self):
        """
        Do not let equally named TIFF/ND2 stems overwrite projections.
        """
        stems = {}
        for name in self.records:
            stems.setdefault(Path(name).stem, []).append(name)
        for names in stems.values():
            if len(names) > 1:
                for name in names:
                    self.records[name].update(
                        Status="FAILED",
                        Stage="Inventory",
                        Error="Duplicate output stem: " + ", ".join(names),
                    )
                _LOG.error("Ambiguous image stems: %s", names)
        self.save()


def active_run():
    return _ACTIVE.get()


@contextmanager
def image_attempt(name, stage, *, final=False):
    run = active_run()
    if run is None:
        yield
    else:
        with run.attempt(name, stage, final=final):
            yield
