#!/usr/bin/env python3
"""Collect the latest valid UMA summaries for each JSON source folder.

Usage: uma_collect_results -i input_paths.json
       python code/collect_results.py -i input_paths.json

Selection is independent for each assay. Image-set mismatches never
trigger a search for older matching runs. CSV bytes are preserved.
Only the Python standard library is used; ImageJ is not started.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from .common.config import read_config
from .common.contracts import (
    ALIGNMENT_SUFFIX,
    ASSAY_TIMESTAMP_PATTERN,
    IMAGE_EXTENSIONS,
    SUMMARY_NAMES,
    THICKNESS_METRICS,
)
from .common.errors import ValidationError
from .common.files import safe_label, save_csv, save_json, sha256_file
from .common.run import close_logger, make_logger, unique_output, utc_now
from .common.version import package_version

SCRIPT_VERSION = "1.0.0"
SELECTION_COLUMNS = [
    "Analysis",
    "Run_Directory",
    "Summary_Path",
    "Run_Timestamp",
    "Status",
    "Reason",
    "Rows",
    "SHA256",
]
CHECK_COLUMNS = [
    "Image_File_Name",
    "Alignment_Present",
    "Alignment_File_Name",
    "Thickness_Present",
    "Thickness_File_Name",
    "Area_Present",
    "Area_File_Name",
    "Check",
    "Details",
]


@dataclass(frozen=True)
class Assay:
    name: str
    prefix: str
    pattern: str
    summary: str
    required: tuple[str, ...]
    numeric: tuple[str, ...]


ASSAYS = (
    Assay(
        "Alignment",
        "Alignment_assay_results_angle_",
        r"Alignment_assay_results_angle_.+_" + ASSAY_TIMESTAMP_PATTERN,
        "Analysis/" + SUMMARY_NAMES["Alignment"],
        (
            "File_Name",
            "Number_of_Z_Stacks",
            "Z_Stack_Type",
            "Orientation_Mode",
        ),
        (),
    ),
    Assay(
        "Thickness",
        "Thickness_assay_results_",
        r"Thickness_assay_results_" + ASSAY_TIMESTAMP_PATTERN,
        SUMMARY_NAMES["Thickness"],
        ("File_Name", *THICKNESS_METRICS),
        THICKNESS_METRICS,
    ),
    Assay(
        "Area",
        "Area_assay_results_",
        r"Area_assay_results_" + ASSAY_TIMESTAMP_PATTERN,
        SUMMARY_NAMES["Area"],
        ("File_Name", "FN_Positive_Pixels", "FN_Area_Percent", "FN_Area"),
        ("FN_Positive_Pixels", "FN_Area_Percent", "FN_Area"),
    ),
)


@dataclass
class Table:
    assay: Assay
    run: Path
    path: Path
    timestamp: datetime
    data: bytes
    rows: list[dict[str, str]]
    status_data: bytes | None

    @property
    def digest(self):
        return hashlib.sha256(self.data).hexdigest()


def create_output(parent, label):
    """Allocate a collection folder with its established UTC name."""
    _, output = unique_output(parent, f"Combined_Results_{label}_")
    return output


def image_key(assay, name):
    """Return a full filename or the stem encoded by alignment."""
    if (
        not name
        or not name.strip()
        or name.startswith(".")
        or "/" in name
        or "\\" in name
        or "\x00" in name
    ):
        raise ValidationError(f"Invalid or hidden image filename: {name!r}")
    if assay.name == "Alignment" and name.endswith(ALIGNMENT_SUFFIX):
        stem = name[: -len(ALIGNMENT_SUFFIX)]
        if not stem or stem.startswith("."):
            raise ValidationError(f"Invalid alignment image stem: {name!r}")
        return "stem", stem
    if Path(name).suffix.lower() not in IMAGE_EXTENSIONS:
        raise ValidationError(
            f"Unrecognized original image filename: {name!r}"
        )
    return "filename", name


def read_summary(assay, run, timestamp):
    path = run / assay.summary
    # Never follow a summary or Analysis directory into a different run.
    if path.is_symlink() or path.parent.is_symlink() or not path.is_file():
        raise ValidationError(
            "Final summary is missing, not a regular file, "
            "or is a symbolic link"
        )
    status_path = run / "run_status.json"
    status_data = None
    status = {}
    if status_path.exists():
        status_data = status_path.read_bytes()
        status = json.loads(status_data.decode("utf-8-sig"))
        if not isinstance(status, dict) or status.get("status") != "SUCCESS":
            state = (
                status.get("status")
                if isinstance(status, dict)
                else "invalid JSON object"
            )
            raise ValidationError(
                f"Run completion status is {state!r}, not SUCCESS"
            )
    data = path.read_bytes()
    reader = csv.DictReader(
        io.StringIO(data.decode("utf-8-sig"), newline=""), strict=True
    )
    header = reader.fieldnames
    if (
        not header
        or len(set(header)) != len(header)
        or any(not key for key in header)
    ):
        raise ValidationError(
            "CSV header is empty or contains duplicate/empty columns"
        )
    missing = sorted(set(assay.required) - set(header))
    if missing:
        raise ValidationError(
            "Required columns are missing: " + ", ".join(missing)
        )
    numeric = list(assay.numeric)
    if assay.name == "Alignment":
        percentages = [
            key
            for key in header
            if re.fullmatch(r"Percentage_Fibers_Aligned_Within_.+_Degree", key)
        ]
        if len(percentages) != 1:
            raise ValidationError(
                "Expected exactly one alignment percentage column"
            )
        numeric.extend(percentages)
    rows, seen = [], set()
    for row in reader:
        line = reader.line_num
        if None in row or any(value is None for value in row.values()):
            raise ValidationError(
                f"CSV row {line} has the wrong number of fields"
            )
        key = image_key(assay, row["File_Name"])
        if key in seen:
            raise ValidationError(
                f"Duplicate image in CSV row {line}: {row['File_Name']}"
            )
        seen.add(key)
        for field in numeric:
            try:
                valid = math.isfinite(float(row[field]))
            except ValueError:
                valid = False
            if not valid:
                raise ValidationError(
                    f"CSV row {line} has an empty or nonfinite {field}"
                )
        if assay.name == "Alignment" and not row["Orientation_Mode"].strip():
            raise ValidationError(
                f"CSV row {line} has an empty Orientation_Mode"
            )
        if (
            assay.name == "Area"
            and "Image_ID" in row
            and row["Image_ID"] != row["File_Name"]
        ):
            raise ValidationError(
                f"CSV row {line}: Image_ID differs from File_Name"
            )
        rows.append(row)
    if not rows:
        raise ValidationError("CSV contains no image rows")
    if assay.name == "Area":
        for field in ("input_images", "processed_images", "generated_masks"):
            if field in status and status[field] != len(rows):
                raise ValidationError(
                    f"Area row count ({len(rows)}) disagrees with "
                    f"run status {field}={status[field]!r}"
                )
    return Table(assay, run, path, timestamp, data, rows, status_data)


def select_latest(source, assay, records, logger):
    """Select by run-name timestamp independently of other assays."""
    candidates = []
    for run in sorted(source.iterdir()):
        if run.name.startswith(".") or not run.name.startswith(assay.prefix):
            continue
        if not run.is_dir() or run.is_symlink():
            continue
        record = dict(zip(SELECTION_COLUMNS, [""] * len(SELECTION_COLUMNS)))
        record.update(
            Analysis=assay.name,
            Run_Directory=str(run),
            Summary_Path=str(run / assay.summary),
            Status="NOT_CHECKED",
        )
        records.append(record)
        try:
            match = re.fullmatch(assay.pattern, run.name)
            if not match:
                raise ValueError(
                    "Run directory does not have a recognized timestamp"
                )
            timestamp = datetime.strptime(match[1], "%Y%m%d_%H%M%S")
            if match.lastindex == 2 and match[2]:
                timestamp = timestamp.replace(microsecond=int(match[2]))
            record["Run_Timestamp"] = timestamp.isoformat()
            candidates.append((timestamp, run, record))
        except ValueError as error:
            record.update(Status="INVALID", Reason=str(error))
            logger.warning("Skipping %s: %s", run, error)
    candidates.sort(key=lambda item: (item[0], item[1].name), reverse=True)
    while candidates:
        latest = candidates[0][0]
        group = [item for item in candidates if item[0] == latest]
        candidates = [item for item in candidates if item[0] != latest]
        valid = []
        for timestamp, run, record in group:
            try:
                table = read_summary(assay, run, timestamp)
                record.update(
                    Status="VALID", Rows=len(table.rows), SHA256=table.digest
                )
                valid.append((table, record))
            except (OSError, ValueError, csv.Error) as error:
                record.update(Status="INVALID", Reason=str(error))
                logger.warning("Skipping %s: %s", run, error)
        if len(valid) > 1:
            for table, record in valid:
                record.update(
                    Status="AMBIGUOUS",
                    Reason="Several valid runs have the same latest timestamp",
                )
            raise ValidationError(
                f"{assay.name}: several valid runs have the same "
                "latest timestamp: "
                + ", ".join(str(table.run) for table, _ in valid)
            )
        if valid:
            table, record = valid[0]
            record.update(
                Status="SELECTED",
                Reason="Latest valid result for this analysis",
            )
            logger.info(
                "Selected %s: %s (%d images)",
                assay.name,
                table.path,
                len(table.rows),
            )
            if table.status_data is None:
                logger.info(
                    "%s has no completion marker; validated its final "
                    "CSV structure and measurements.",
                    table.run,
                )
            for _, _, older in candidates:
                older["Reason"] = "Older than the selected valid result"
            return table
    return None


def compare_images(source, tables):
    """Resolve alignment stems, retaining other assays' extensions."""
    originals = {
        path.name
        for path in source.iterdir()
        if not path.name.startswith(".")
        and path.is_file()
        and path.suffix.lower() in IMAGE_EXTENSIONS
    }
    for table in tables.values():
        for row in table.rows:
            kind, name = image_key(table.assay, row["File_Name"])
            if kind == "filename":
                originals.add(name)
    by_stem = {}
    for name in originals:
        by_stem.setdefault(Path(name).stem, set()).add(name)
    maps, errors, unresolved = {}, [], []
    for analysis, table in tables.items():
        mapping = {}
        maps[analysis] = mapping
        for row in table.rows:
            source_name = row["File_Name"]
            kind, name = image_key(table.assay, source_name)
            if kind == "stem":
                choices = sorted(by_stem.get(name, set()))
                if len(choices) != 1:
                    detail = (
                        f"Alignment stem {name!r} is ambiguous: {choices}"
                        if choices
                        else (
                            "Cannot recover the original filename for "
                            f"alignment stem {name!r}"
                        )
                    )
                    errors.append(detail)
                    unresolved.append(
                        {
                            "Image_File_Name": "",
                            "Alignment_Present": True,
                            "Alignment_File_Name": source_name,
                            "Check": "UNRESOLVED",
                            "Details": detail,
                        }
                    )
                    continue
                name = choices[0]
            if name in mapping:
                errors.append(f"{analysis}: duplicate resolved image {name!r}")
            mapping[name] = source_name
    all_images = sorted(
        set().union(*(set(mapping) for mapping in maps.values()))
    )
    report = []
    for name in all_images:
        absent = [
            analysis
            for analysis, mapping in maps.items()
            if name not in mapping
        ]
        check = (
            "MISMATCH"
            if absent
            else ("MATCH" if len(tables) > 1 else "NOT_COMPARABLE")
        )
        detail = "Missing from selected " + ", ".join(absent) if absent else ""
        if absent:
            errors.append(f"{name}: {detail}")
        record = {"Image_File_Name": name, "Check": check, "Details": detail}
        for assay in ASSAYS:
            mapping = maps.get(assay.name)
            record[assay.name + "_Present"] = (
                name in mapping if mapping is not None else "NOT_AVAILABLE"
            )
            record[assay.name + "_File_Name"] = (
                mapping.get(name, "") if mapping is not None else ""
            )
        report.append(record)
    report.extend(unresolved)
    return report, errors


def remove_copies(output, label):
    for assay in ASSAYS:
        target = output / f"{label}_{Path(assay.summary).name}"
        target.unlink(missing_ok=True)
        target.with_name(target.name + ".partial").unlink(missing_ok=True)


def publish_copies(output, label, tables):
    """Publish validated snapshots while their sources remain stable."""
    copies = []
    try:
        for table in tables.values():
            target = output / f"{label}_{table.path.name}"
            temporary = target.with_name(target.name + ".partial")
            temporary.write_bytes(table.data)
            if sha256_file(temporary) != table.digest:
                raise OSError(f"Copy verification failed: {temporary}")
            copies.append(
                {
                    "analysis": table.assay.name,
                    "path": str(target),
                    "sha256": table.digest,
                }
            )
        for table in tables.values():
            status_path = table.run / "run_status.json"
            current_status = (
                status_path.read_bytes() if status_path.exists() else None
            )
            if (
                table.path.read_bytes() != table.data
                or current_status != table.status_data
            ):
                raise ValidationError(
                    f"Selected source changed during collection: {table.path}"
                )
        for item in copies:
            target = Path(item["path"])
            target.with_name(target.name + ".partial").replace(target)
        return copies
    except BaseException:
        remove_copies(output, label)
        raise


def collect_folder(source, input_json):
    label = safe_label(source.name)
    # Do not create missing or unwritable input folders.
    try:
        output = create_output(source, label)
    except OSError as error:
        return diagnostic_failure(source, input_json, error)
    logger = make_logger(output)
    status = {
        "status": "RUNNING",
        "script_version": SCRIPT_VERSION,
        "package_version": package_version(),
        "started_utc": utc_now(),
        "input_json": str(input_json),
        "source_folder": str(source),
        "source_name": source.name,
        "output_label": label,
        "run_directory": str(output),
        "selection_rule": (
            "Latest valid result per assay, independently; run-name timestamp"
        ),
        "copied_csvs": [],
        "errors": [],
    }
    records, tables, selection_errors = [], {}, []
    try:
        save_json(output / "run_status.json", status)
        logger.info(
            "UMA results collector %s (package %s). Source: %s",
            SCRIPT_VERSION,
            package_version(),
            source,
        )
        logger.info("Output: %s", output)
        for assay in ASSAYS:
            try:
                table = select_latest(source, assay, records, logger)
                if table is not None:
                    tables[assay.name] = table
            except ValidationError as error:
                selection_errors.append(str(error))
        save_csv(output / "selection_report.csv", SELECTION_COLUMNS, records)
        missing = [assay.name for assay in ASSAYS if assay.name not in tables]
        status.update(
            analyses_found=len(tables),
            analyses_expected=len(ASSAYS),
            missing_analyses=missing,
            selected={
                name: {
                    "summary": str(table.path),
                    "run_directory": str(table.run),
                    "run_timestamp": table.timestamp.isoformat(),
                    "sha256": table.digest,
                    "image_rows": len(table.rows),
                }
                for name, table in tables.items()
            },
        )
        message = f"Found {len(tables)} of 3 analyses."
        if missing:
            logger.warning(
                "%s Missing or without a selectable valid result: %s",
                message,
                ", ".join(missing),
            )
        else:
            logger.info(message)
        report, errors = compare_images(source, tables)
        errors = selection_errors + errors
        save_csv(output / "image_check.csv", CHECK_COLUMNS, report)
        if not tables:
            errors.append("No valid analysis summaries were found")
        if errors:
            for error in errors:
                logger.error(error)
            status.update(
                status="VALIDATION_FAILED", image_check="FAILED", errors=errors
            )
            logger.error(
                "Validation failed. Only diagnostics were saved; "
                "no summary CSVs were collected."
            )
        else:
            status["image_check"] = (
                "MATCH" if len(tables) > 1 else "NOT_COMPARABLE"
            )
            if len(tables) == 1:
                logger.warning(
                    "Only one analysis is available; comparison between "
                    "analyses is not possible."
                )
            else:
                logger.info(
                    "Image sets match across all %d selected analyses "
                    "(%d images).",
                    len(tables),
                    len(report),
                )
            status["copied_csvs"] = publish_copies(output, label, tables)
            status["status"] = (
                "SUCCESS_WITH_MISSING_ANALYSES" if missing else "SUCCESS"
            )
            logger.info(
                "Collected %d unchanged summary CSV(s). Output: %s",
                len(tables),
                output,
            )
        status["ended_utc"] = utc_now()
        save_json(output / "run_status.json", status)
        return status["status"].startswith("SUCCESS"), output
    except (Exception, KeyboardInterrupt) as error:
        remove_copies(output, label)
        status.update(
            status="CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else "ERROR",
            ended_utc=utc_now(),
            errors=[str(error) or "Interrupted"],
            copied_csvs=[],
        )
        logger.exception("Collection failed; no summary CSVs were collected.")
        (output / "traceback.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        save_csv(output / "selection_report.csv", SELECTION_COLUMNS, records)
        save_json(output / "run_status.json", status)
        if isinstance(error, KeyboardInterrupt):
            raise
        return False, output
    finally:
        close_logger(logger)


def diagnostic_failure(source, input_json, error):
    """Keep CWD startup diagnostics if source output is impossible."""
    output = None
    try:
        label = (
            safe_label(source.name)
            if source is not None
            else "configuration_error"
        )
        output = create_output(Path.cwd(), label)
        logger = make_logger(output)
        try:
            logger.error("Collection could not start: %s", error)
            logger.info("Diagnostics: %s", output)
            save_json(
                output / "run_status.json",
                {
                    "status": "ERROR",
                    "script_version": SCRIPT_VERSION,
                    "ended_utc": utc_now(),
                    "source_folder": str(source)
                    if source is not None
                    else None,
                    "input_json": str(input_json),
                    "errors": [str(error)],
                    "copied_csvs": [],
                },
            )
        finally:
            close_logger(logger)
    except OSError as diagnostic_error:
        print(
            f"ERROR: {error}. Could not save diagnostics: {diagnostic_error}",
            file=sys.stderr,
        )
    return False, output


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Collect and verify UMA alignment, thickness, and area summaries"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="JSON file containing folder_paths",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {package_version()} (collector {SCRIPT_VERSION})",
    )
    args = parser.parse_args(argv)
    input_json = Path(args.input).expanduser().absolute()
    try:
        folders = read_config(input_json)
    except (OSError, ValueError) as error:
        diagnostic_failure(None, input_json, error)
        return 1
    seen, failed, succeeded = set(), 0, 0
    for source in folders:
        try:
            canonical = source.resolve()
        except (OSError, ValueError) as error:
            diagnostic_failure(source, input_json, error)
            failed += 1
            continue
        if canonical in seen:
            print(f"Skipping duplicate JSON folder: {source}", flush=True)
            continue
        seen.add(canonical)
        try:
            success, _ = collect_folder(source, input_json)
        except KeyboardInterrupt:
            print("Collection interrupted.", file=sys.stderr)
            return 130
        except Exception as error:
            diagnostic_failure(source, input_json, error)
            success = False
        succeeded += int(success)
        failed += int(not success)
    print(
        f"Collection finished: {succeeded} folder(s) succeeded; "
        f"{failed} failed.",
        flush=True,
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
