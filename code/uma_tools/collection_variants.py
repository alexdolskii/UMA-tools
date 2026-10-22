"""Collect independent Area variants and create immutable revisions."""

from __future__ import annotations

import csv
import hashlib
import re
import sys
import traceback
from datetime import datetime
from pathlib import Path

from . import collect_results as core
from . import package_version, report_inputs
from .area_variants import (
    COLLECTION_SCHEMA,
    inside,
    profile,
    sort_key,
    variant_id,
    variant_label,
)
from .config import ValidationError
from .files import (
    assay_directory,
    safe_label,
    save_csv,
    save_json,
    sha256_file,
)
from .progress import console, folder_logged, outcome, phase
from .report_options import _ask, _choice
from .report_schema import ValidationError as ReportValidationError
from .run import close_logger, make_logger, utc_now


class SelectionLog:
    def __init__(self, logger):
        self.logger = logger

    def event(self, level, stage, message, **kwargs):
        self.logger.info("%s: %s: %s", level, stage, message)


class DiscoveryLog:
    """Keep discovery diagnostics until a revision directory is allocated."""

    def __init__(self):
        self.messages = []

    def info(self, message, *args):
        self.messages.append(("info", message % args if args else message))

    def warning(self, message, *args):
        self.messages.append(("warning", message % args if args else message))

    def replay(self, logger):
        for level, message in self.messages:
            getattr(logger, level)(message)


def verify_tables(tables):
    for table in tables:
        path = table.run / "run_status.json"
        state = path.read_bytes() if path.exists() else None
        if table.path.read_bytes() != table.data or state != table.status_data:
            raise ValidationError(
                f"Selected source changed during collection: {table.path}"
            )


def previous_collection(source, logger):
    if not assay_directory(source).is_dir():
        return None
    try:
        return report_inputs.select_combined(source, SelectionLog(logger))
    except ReportValidationError as error:
        if str(error).startswith("No successful Combined_Results"):
            return None
        raise


def scan_area(source, records, logger):
    """Deduplicate only within an exact profile; reject timestamp ties."""
    groups = {}
    assay = core.ASSAYS[2]
    for run in sorted(assay_directory(source).iterdir()):
        if (
            not run.name.startswith(assay.prefix)
            or run.is_symlink()
            or not run.is_dir()
        ):
            continue
        record = dict.fromkeys(core.SELECTION_COLUMNS, "")
        record.update(
            Analysis="Area",
            Run_Directory=str(run),
            Summary_Path=str(run / assay.summary),
            Status="INVALID",
        )
        records.append(record)
        try:
            match = re.fullmatch(assay.pattern, run.name)
            if match is None:
                raise ValidationError("Unrecognized Area run timestamp")
            timestamp = datetime.strptime(match[1], "%Y%m%d_%H%M%S")
            if match[2]:
                timestamp = timestamp.replace(microsecond=int(match[2]))
            table = core.read_summary(assay, run, timestamp)
            record.update(Rows=len(table.rows), SHA256=table.digest)
            settings = profile(table.rows)
            key = variant_id(settings)
            record.update(
                Run_Timestamp=timestamp.isoformat(),
                Status="VALID",
                Rows=len(table.rows),
                SHA256=table.digest,
            )
            groups.setdefault(key, []).append((table, settings, record))
        except (OSError, ValueError, csv.Error) as error:
            record["Reason"] = str(error)
            logger.warning("Skipping %s: %s", run, error)
    chosen, errors = {}, []
    for key, candidates in groups.items():
        candidates.sort(key=lambda item: item[0].timestamp, reverse=True)
        latest = candidates[0][0].timestamp
        tied = [item for item in candidates if item[0].timestamp == latest]
        if len(tied) > 1:
            message = (
                f"Area {variant_label(tied[0][1])}: several valid runs have "
                "the same latest timestamp"
            )
            errors.append(message)
            for _, _, record in tied:
                record.update(Status="AMBIGUOUS", Reason=message)
            continue
        table, settings, record = tied[0]
        record.update(Status="SELECTED", Reason="Latest valid run for profile")
        for _, _, older in candidates[1:]:
            older["Reason"] = "Older run of the same Area profile"
        chosen[key] = (table, settings, None)
    return chosen, errors


def archive_table(table, output, relative, origin=None):
    """Preserve complete rows and their audit before variant exclusions."""
    folder = output / relative
    path = folder / table.assay.summary
    path.parent.mkdir(parents=True, exist_ok=True)
    if table.path.read_bytes() != table.data:
        raise ValidationError(
            f"Source changed during collection: {table.path}"
        )
    status_path = table.run / "run_status.json"
    status_data = status_path.read_bytes() if status_path.exists() else None
    if status_data != table.status_data:
        raise ValidationError(f"Source status changed: {status_path}")
    path.write_bytes(table.data)
    if sha256_file(path) != table.digest:
        raise OSError(f"Copy verification failed: {path}")
    if status_data is not None:
        (folder / "run_status.json").write_bytes(status_data)
    origin = origin or {}
    record = {
        "analysis": table.assay.name,
        "snapshot": relative,
        "sha256": table.digest,
        "status_sha256": (
            hashlib.sha256(status_data).hexdigest() if status_data else None
        ),
        "timestamp": table.timestamp.isoformat(),
        "origin_run": origin.get("origin_run", str(table.run)),
        "origin_summary": origin.get("origin_summary", str(table.path)),
    }
    return load_snapshot(output, record), record


def load_snapshot(collection, record):
    assay = next(
        (item for item in core.ASSAYS if item.name == record.get("analysis")),
        None,
    )
    if assay is None:
        raise ValidationError("Unknown frozen assay")
    folder = inside(collection, record["snapshot"])
    table = core.read_summary(
        assay, folder, datetime.fromisoformat(record["timestamp"])
    )
    status_hash = (
        hashlib.sha256(table.status_data).hexdigest()
        if table.status_data
        else None
    )
    if (
        table.digest != record["sha256"]
        or status_hash != record["status_sha256"]
    ):
        raise ValidationError(f"Frozen input changed: {folder}")
    return table


def selection_record(table, origin):
    return {
        "summary": origin["origin_summary"],
        "run_directory": origin["origin_run"],
        "run_timestamp": table.timestamp.isoformat(),
        "sha256": table.digest,
        "image_rows": len(table.rows),
    }


def build_variant(source, output, tables, origins, settings):
    """Keep the existing single-report contract inside each variant."""
    output.mkdir(parents=True, exist_ok=True)
    report, errors = core.compare_images(source, tables)
    save_csv(output / "image_check.csv", core.CHECK_COLUMNS, report)
    exclusions = [
        {"Image_File_Name": row["Image_File_Name"], "Reason": row["Details"]}
        for row in report
        if row["Check"] == "EXCLUDED"
    ]
    retained = [row for row in report if row["Check"] == "MATCH"]
    save_csv(
        output / "processing_exclusions.csv",
        ["Image_File_Name", "Reason"],
        exclusions,
    )
    missing = [assay.name for assay in core.ASSAYS if assay.name not in tables]
    if not retained and not missing:
        errors.append("No matched successful images remain")
    state = {
        "status": "RUNNING",
        "source_folder": str(source),
        "source_name": source.name,
        "copied_csvs": [],
        "analyses_found": len(tables),
        "analyses_expected": 3,
        "missing_analyses": missing,
        "area_profile": settings,
        "selected": {
            name: selection_record(table, origins[name])
            for name, table in tables.items()
        },
        "image_check": "FAILED" if errors else "MATCH",
        "excluded_images": len(exclusions),
        "retained_images": len(retained),
        "processing_exclusions": exclusions,
        "exclusions_sha256": sha256_file(output / "processing_exclusions.csv"),
        "errors": errors,
    }
    if errors:
        state["status"] = "VALIDATION_FAILED"
    else:
        state["copied_csvs"] = core.publish_copies(
            output, safe_label(source.name), tables, image_report=report
        )
        state["status"] = (
            "SUCCESS_WITH_MISSING_ANALYSES" if missing else "SUCCESS"
        )
    state["ended_utc"] = utc_now()
    save_json(output / "run_status.json", state)
    return state


def copy_template(previous, output):
    if previous is None:
        return None
    templates = [
        path
        for path in previous.iterdir()
        if path.suffix.lower() == ".xlsx"
        and not path.name.startswith((".", "~$"))
        and path.is_file()
        and not path.is_symlink()
    ]
    if len(templates) > 1:
        raise ValidationError(
            "Several plate templates in the previous collection: "
            + ", ".join(path.name for path in templates)
        )
    if not templates:
        return None
    path = templates[0]
    content = path.read_bytes()
    target = output / path.name
    target.write_bytes(content)
    if target.read_bytes() != content or path.read_bytes() != content:
        raise ValidationError("Plate template changed during collection")
    return {"file": target.name, "sha256": sha256_file(target)}


def plan(source, previous, previous_state, records, logger, args):
    areas, errors = scan_area(source, records, logger)
    modern = previous_state.get("collection_schema") == COLLECTION_SCHEMA
    old_areas = previous_state.get("area_runs", {}) if modern else {}
    changes = []
    for key, (table, settings, _) in areas.items():
        old = old_areas.get(key)
        if old is None:
            changes.append(f"NEW: {variant_label(settings)}")
        elif table.timestamp < datetime.fromisoformat(old["timestamp"]):
            # Removing a source run must not downgrade a frozen revision.
            continue
        elif (
            old["origin_run"] != str(table.run)
            or old["sha256"] != table.digest
            or old["status_sha256"]
            != (
                hashlib.sha256(table.status_data).hexdigest()
                if table.status_data
                else None
            )
        ):
            changes.append(f"UPDATED: {variant_label(settings)}")
    if modern:
        for variant in previous_state.get("variants", []):
            if variant.get("status") == "ERROR":
                changes.append(f"RETRY: {variant['label']}")
    console(f"{source.name}: {len(areas)} valid Area profile(s).")
    for change in changes:
        console(change)
    mode = getattr(args, "collection_mode", None)
    interactive = not getattr(args, "non_interactive", False)
    if (
        mode is None
        and previous is not None
        and sys.stdin.isatty()
        and interactive
    ):
        console(f"Previous collection: {previous}")
        mode = _ask(
            "Collection: update Area only, rebuild all, or keep",
            "update" if modern else "rebuild",
            _choice({key: key for key in ("update", "rebuild", "keep")}),
        )
    mode = mode or ("update" if modern else "rebuild")
    if mode == "keep":
        if previous is None:
            raise ValidationError("No completed collection to keep")
        return None
    if mode == "update" and not modern:
        console("First multi-Area collection: rebuilding saved assay tables.")
        mode = "rebuild"
    shared, origins = {}, {}
    if mode == "update":
        for name, record in previous_state["shared_inputs"].items():
            shared[name] = load_snapshot(previous, record)
            origins[name] = record
        for key, record in old_areas.items():
            old_table = load_snapshot(previous, record)
            settings = profile(old_table.rows)
            if variant_id(settings) != key:
                raise ValidationError("Frozen Area profile ID disagrees")
            new = areas.get(key)
            if new is None or new[0].timestamp < old_table.timestamp:
                areas[key] = (old_table, settings, record)
            elif (
                new[0].timestamp == old_table.timestamp
                and str(new[0].run) != record["origin_run"]
            ):
                raise ValidationError("Ambiguous Area revision timestamp")
        if not changes and not errors:
            outcome(
                previous_state["status"],
                f"{source}: no new Area; collection unchanged.",
            )
            return None
    else:
        for assay in core.ASSAYS[:2]:
            table = core.select_latest(source, assay, records, logger)
            if table is not None:
                shared[assay.name] = table
    return shared, origins, areas, errors, mode, changes


@folder_logged("source")
def collect(source, input_json, args):
    """Publish only completed revisions; never modify a previous folder."""
    records = []
    discovery = DiscoveryLog()
    logger = discovery
    output = None
    state = {}
    try:
        if not source.is_dir():
            raise ValidationError(f"Source directory does not exist: {source}")
        assay_directory(source, create=True)
        previous = previous_collection(source, logger)
        previous_state = (
            report_inputs.read_collection_status(previous) if previous else {}
        )
        planned = plan(source, previous, previous_state, records, logger, args)
        if planned is None:
            return core.CollectStatus(previous_state.get("status", "SUCCESS"))
        shared, origins, areas, errors, mode, changes = planned
        if not areas:
            if any(
                record["Rows"] and record["Status"] == "INVALID"
                for record in records
            ):
                raise ValidationError(
                    "No usable Area profile; check recorded threshold settings"
                )
            # Preserve useful missing-assay diagnostics for legacy datasets.
            return core.collect_folder(source, input_json)
        output = core.create_output(
            assay_directory(source), safe_label(source.name)
        )
        logger = make_logger(output)
        discovery.replay(logger)
        watched = list(shared.values()) + [item[0] for item in areas.values()]
        state = {
            "collection_schema": COLLECTION_SCHEMA,
            "status": "RUNNING",
            "package_version": package_version(),
            "source_folder": str(source),
            "source_name": source.name,
            "run_directory": str(output),
            "input_json": str(input_json),
            "started_utc": utc_now(),
            "mode": mode,
            "parent_collection": str(previous) if previous else None,
            "selection_rule": (
                "Latest valid run per Area profile; "
                "frozen shared assays on update"
            ),
            "changes": changes,
            "errors": errors,
            "shared_inputs": {},
            "area_runs": {},
            "variants": [],
            "copied_csvs": [],
        }
        save_json(output / "run_status.json", state)
        logger.info("Collection mode: %s; parent: %s", mode, previous)
        for change in changes:
            logger.info(change)
        for name, table in list(shared.items()):
            shared[name], state["shared_inputs"][name] = archive_table(
                table, output, f"Raw/{name}", origins.get(name)
            )
        state["template"] = copy_template(previous, output)
        for index, (key, (table, settings, origin)) in enumerate(
            sorted(areas.items(), key=lambda item: sort_key(item[1][1])), 1
        ):
            label = variant_label(settings)
            phase(f"Area variant {index}/{len(areas)}: {label}")
            item = {
                "id": key,
                "label": label,
                "profile": settings,
                "path": f"Variants/{label}_{key}",
                "status": "RUNNING",
            }
            state["variants"].append(item)
            try:
                frozen, area_origin = archive_table(
                    table, output, f"Raw/Area/{key}", origin
                )
                state["area_runs"][key] = area_origin
                variant = build_variant(
                    source,
                    output / item["path"],
                    {**shared, "Area": frozen},
                    {**state["shared_inputs"], "Area": area_origin},
                    settings,
                )
                item.update(
                    status=variant["status"],
                    errors=variant["errors"],
                    excluded_images=variant["excluded_images"],
                    retained_images=variant["retained_images"],
                )
                logger.info(
                    "%s: %s; %s", label, item["status"], item["errors"]
                )
            except (OSError, ValueError, csv.Error) as error:
                item.update(status="ERROR", errors=[str(error)])
                logger.exception("Variant failed: %s", label)
                diagnostic = output / item["path"]
                diagnostic.mkdir(parents=True, exist_ok=True)
                trace = traceback.format_exc()
                (diagnostic / "traceback.txt").write_text(trace)
                with (output / "traceback.txt").open("a") as stream:
                    stream.write(f"{label}\n{trace}\n")
        good = [
            v for v in state["variants"] if v["status"].startswith("SUCCESS")
        ]
        failed = len(state["variants"]) - len(good)
        state["status"] = (
            "VALIDATION_FAILED"
            if not good
            else "PARTIAL"
            if failed or errors
            else "SUCCESS_WITH_MISSING_ANALYSES"
            if len(shared) < 2
            else "SUCCESS"
        )
        # Keep the familiar flat CSVs for a single-profile collection.
        # Complete, unfiltered inputs always remain in Raw/.
        verify_tables(watched)
        only = state["variants"][0] if len(state["variants"]) == 1 else None
        if only and (output / only["path"] / "run_status.json").is_file():
            child = output / only["path"]
            single = report_inputs.read_collection_status(child)
            state["errors"].extend(single["errors"])
            for filename in ("image_check.csv", "processing_exclusions.csv"):
                (output / filename).write_bytes(
                    (child / filename).read_bytes()
                )
            for entry in single["copied_csvs"]:
                original = Path(entry["path"])
                target = output / original.name
                target.write_bytes(original.read_bytes())
                entry["path"] = str(target)
            for name in (
                "copied_csvs",
                "selected",
                "analyses_found",
                "analyses_expected",
                "missing_analyses",
                "image_check",
                "excluded_images",
                "retained_images",
                "processing_exclusions",
                "exclusions_sha256",
            ):
                state[name] = single[name]
        state["ended_utc"] = utc_now()
        save_csv(
            output / "selection_report.csv", core.SELECTION_COLUMNS, records
        )
        save_json(output / "run_status.json", state)
        outcome(
            state["status"],
            f"{source}: {len(shared) + 1} analyses; "
            f"{len(good)} usable Area variant(s); "
            f"{failed} failed. Output: {output}",
        )
        return core.CollectStatus(state["status"])
    except (Exception, KeyboardInterrupt, EOFError) as error:
        if output is None:
            if isinstance(error, (KeyboardInterrupt, EOFError)):
                raise KeyboardInterrupt from error
            try:
                output = core.create_output(
                    assay_directory(source), safe_label(source.name)
                )
                logger = make_logger(output)
                discovery.replay(logger)
                state = {
                    "source_folder": str(source),
                    "copied_csvs": [],
                    "collection_schema": COLLECTION_SCHEMA,
                }
            except OSError:
                return core.diagnostic_failure(source, input_json, error)
        state.update(
            status="CANCELLED"
            if isinstance(error, (KeyboardInterrupt, EOFError))
            else "ERROR",
            errors=[str(error)],
            ended_utc=utc_now(),
        )
        save_json(output / "run_status.json", state)
        (output / "traceback.txt").write_text(traceback.format_exc())
        logger.exception("Collection revision failed; previous data preserved")
        if state["status"] == "CANCELLED":
            raise KeyboardInterrupt from error
        return core.CollectStatus.ERROR
    finally:
        if output is not None:
            close_logger(logger)
