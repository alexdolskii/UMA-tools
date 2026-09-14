"""
Select and archive collected inputs, with report file and log policies.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime
from functools import partial
from pathlib import Path, PureWindowsPath

from .contracts import SUMMARY_NAMES
from .files import save_csv as _save_csv
from .files import save_json as _save_json
from .files import sha256_file
from .report_schema import EventLogger, ReportInputs, ValidationError
from .run import RunLog as RunLog
from .run import utc_now as _utc_now

save_csv = partial(_save_csv, encoding="utf-8-sig")
save_json = partial(
    _save_json, atomic=False, allow_nan=False, trailing_newline=False
)
utc_now = partial(_utc_now, timespec="seconds")

ROLES = {
    "Alignment": ("alignment", SUMMARY_NAMES["Alignment"]),
    "Thickness": ("thickness", SUMMARY_NAMES["Thickness"]),
    "Area": ("fibronectin", SUMMARY_NAMES["Area"]),
}
COMBINED_PATTERN = re.compile(
    r"Combined_Results_.+_(\d{8}_\d{6})(?:_(\d{6}))?(?:_\d+)?"
)
SUCCESS_STATES = {"SUCCESS", "SUCCESS_WITH_MISSING_ANALYSES"}


def read_collection_status(combined):
    path = combined / "run_status.json"
    if path.is_symlink() or not path.is_file():
        raise ValidationError(
            "Collector completion record is missing or not a regular file: "
            f"{path}"
        )
    status = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(status, dict):
        raise ValidationError(f"Invalid collector completion record: {path}")
    return status


def select_combined(source: Path, log: EventLogger) -> Path:
    """
    Choose the newest completed collection, without inspecting its
    template.
    """
    candidates = []
    for path in sorted(source.iterdir()):
        if path.name.startswith(".") or not path.name.startswith(
            "Combined_Results_"
        ):
            continue
        if path.is_symlink() or not path.is_dir():
            continue
        try:
            match = COMBINED_PATTERN.fullmatch(path.name)
            if not match:
                raise ValidationError(
                    "Unrecognized Combined_Results timestamp"
                )
            timestamp = datetime.strptime(match[1], "%Y%m%d_%H%M%S")
            if match[2]:
                timestamp = timestamp.replace(microsecond=int(match[2]))
            status = read_collection_status(path)
            if status.get("status") not in SUCCESS_STATES:
                raise ValidationError(
                    f"Collector status is {status.get('status')!r}"
                )
            candidates.append((timestamp, path))
        except (OSError, ValueError, ValidationError) as error:
            log.event(
                "WARNING", "Collection selection", f"Skipping {path}: {error}"
            )
    if not candidates:
        raise ValidationError(
            f"No successful Combined_Results directory found in {source}"
        )
    candidates.sort(key=lambda item: (item[0], item[1].name), reverse=True)
    latest = [
        path for timestamp, path in candidates if timestamp == candidates[0][0]
    ]
    if len(latest) != 1:
        raise ValidationError(
            "Several successful collections have the same latest timestamp: "
            + ", ".join(str(path) for path in latest)
        )
    log.event(
        "INFO",
        "Collection selection",
        f"Selected latest successful collection: {latest[0]}",
    )
    return latest[0]


def recorded_basename(value):
    """
    Use archived filenames even after a collection is moved between
    drives.
    """
    if not isinstance(value, str) or not value:
        raise ValidationError("Collector manifest contains an empty CSV path")
    name = PureWindowsPath(value).name if "\\" in value else Path(value).name
    if not name or name.startswith((".", "~$")) or name in {".", ".."}:
        raise ValidationError(f"Invalid collected CSV filename: {value!r}")
    return name


def collection_entries(status):
    """
    Validate the complete collector manifest, including archived copies.
    """
    if not isinstance(status, dict):
        raise ValidationError("Invalid collector completion record")
    if status.get("status") not in SUCCESS_STATES:
        raise ValidationError(
            "The selected collection is no longer marked successful"
        )
    entries = status.get("copied_csvs")
    if not isinstance(entries, list):
        raise ValidationError(
            "Collector completion record has no copied_csvs manifest"
        )
    by_analysis = {}
    for entry in entries:
        if not isinstance(entry, dict) or entry.get("analysis") not in ROLES:
            raise ValidationError(
                "Collector manifest contains an unknown analysis entry"
            )
        analysis = entry["analysis"]
        if analysis in by_analysis:
            raise ValidationError(
                f"Duplicate {analysis} entries in collector manifest"
            )
        by_analysis[analysis] = entry
    missing = [analysis for analysis in ROLES if analysis not in by_analysis]
    if missing:
        raise ValidationError(
            "All three analyses are required. Missing: " + ", ".join(missing)
        )
    return by_analysis


def discover_inputs(combined: Path, source_name: str) -> ReportInputs:
    """
    Require all three collected summaries and exactly one visible plate
    map.
    """
    by_analysis = collection_entries(read_collection_status(combined))
    paths, details = {}, []
    for analysis, (role, suffix) in ROLES.items():
        entry = by_analysis[analysis]
        name = recorded_basename(entry.get("path"))
        path = combined / name
        digest = entry.get("sha256")
        issue = None
        if not (name == suffix or name.endswith("_" + suffix)):
            issue = f"Unexpected {analysis} summary filename: {name}"
        elif path.is_symlink() or not path.is_file():
            issue = (
                f"Collected {analysis} CSV is missing or not a regular file: "
                f"{path}"
            )
        elif not isinstance(digest, str) or not re.fullmatch(
            r"[0-9a-fA-F]{64}", digest
        ):
            issue = (
                f"Missing or invalid {analysis} SHA256 in collector manifest"
            )
        elif sha256_file(path) != digest.lower():
            issue = (
                f"Collected {analysis} CSV has changed since collection: "
                f"{path}"
            )
        if issue:
            details.append({"Input": role, "Path": str(path), "Issue": issue})
        else:
            paths[role] = path
    templates = sorted(
        path
        for path in combined.iterdir()
        if path.suffix.lower() == ".xlsx"
        and not path.name.startswith((".", "~$"))
        and path.is_file()
        and not path.is_symlink()
    )
    if len(templates) != 1:
        issue = (
            (
                "No plate-template .xlsx found in the selected "
                "Combined_Results directory"
            )
            if not templates
            else "More than one plate-template .xlsx found: "
            + ", ".join(path.name for path in templates)
        )
        details.append(
            {"Input": "template", "Path": str(combined), "Issue": issue}
        )
    else:
        paths["template"] = templates[0]
    if details:
        raise ValidationError(
            "Selected collection is not ready for reporting. "
            + "; ".join(item["Issue"] for item in details),
            details,
        )
    return paths


def archive_inputs(paths, input_json, combined, directory):
    """
    Make verified input snapshots before parsing or plotting their
    contents.
    """
    archive = directory / "Inputs"
    archive.mkdir()
    sources = {
        **paths,
        "configuration": input_json,
        "collector_status": combined / "run_status.json",
    }
    for role, name in (
        ("collector_image_check", "image_check.csv"),
        ("collector_selection_report", "selection_report.csv"),
    ):
        path = combined / name
        if path.is_file() and not path.is_symlink():
            sources[role] = path
    names = {
        "configuration": "input_paths.json",
        "collector_status": "collector_run_status.json",
        "collector_image_check": "collector_image_check.csv",
        "collector_selection_report": "collector_selection_report.csv",
    }
    records, snapshots, used = [], {}, set()
    for role, source in sources.items():
        name = names.get(role, source.name)
        if name in used:
            raise ValidationError(f"Input archive filename collision: {name}")
        used.add(name)
        data = source.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        target = archive / name
        target.write_bytes(data)
        if sha256_file(target) != digest:
            raise OSError(f"Input copy verification failed: {target}")
        records.append(
            {
                "Input": role,
                "Path": str(source),
                "Archived_Path": str(target),
                "Archive_Relative_Path": str(target.relative_to(directory)),
                "Bytes": len(data),
                "SHA256": digest,
            }
        )
        if role in paths:
            snapshots[role] = target
    save_csv(directory / "input_manifest.csv", list(records[0]), records)
    # Reconcile copied summary bytes against the archived collector
    # record as well, closing the gap between discovery and archiving.
    archived_status = json.loads(
        (archive / "collector_run_status.json").read_text(encoding="utf-8-sig")
    )
    entries = collection_entries(archived_status)
    archived_records = {record["Input"]: record for record in records}
    for analysis, (role, _) in ROLES.items():
        entry, record = entries[analysis], archived_records[role]
        digest = entry.get("sha256")
        if (
            recorded_basename(entry.get("path")) != paths[role].name
            or not isinstance(digest, str)
            or digest.lower() != record["SHA256"]
        ):
            raise ValidationError(
                "Collector manifest changed while inputs were copied: "
                + record["Path"]
            )
    verify_sources(records)
    return snapshots, records


def verify_sources(records):
    for record in records:
        if sha256_file(Path(record["Path"])) != record["SHA256"]:
            raise ValidationError(
                "An input changed during report generation: " + record["Path"]
            )


def save_details(path, rows):
    """
    Write every diagnostic field without dropping affected records.
    """
    columns = list(dict.fromkeys(key for row in rows for key in row))
    save_csv(path, columns or ["Issue"], rows)
