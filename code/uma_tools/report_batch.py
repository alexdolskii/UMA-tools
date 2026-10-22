"""Plan independent reports and reuse only verified, unchanged outputs."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json

from . import package_version, report_inputs
from .area_variants import (
    COLLECTION_SCHEMA,
    digest,
    inside,
    profile,
    variant_id,
    variant_label,
)
from .config import ValidationError
from .files import assay_directory, save_csv, save_json, sha256_file
from .progress import console, folder_logged, outcome, phase
from .report_options import _ask, _choice
from .report_schema import ValidationError as ReportValidationError

PLAN_COLUMNS = [
    "Variant_ID",
    "Area_variant",
    "Action",
    "Previous_report",
    "Reason",
]


def identity(source, args, paths, variant_state, variant, fia):
    versions = {}
    for name in ("numpy", "matplotlib", "openpyxl", "Pillow", "scipy"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return {
        "schema": 1,
        "package_version": package_version(),
        "dependencies": versions,
        "source": str(source.resolve()),
        "source_name": source.name,
        "area_profile": variant["profile"],
        "inputs": {key: sha256_file(path) for key, path in paths.items()},
        "selected_runs": variant_state.get("selected", {}),
        "exclusions": variant_state.get("processing_exclusions", []),
        "settings": {
            key: getattr(args, key, default)
            for key, default in (
                ("fn_threshold", 20),
                ("min_nuclei", 0),
                ("with_fia_nuclei", False),
                ("stats_unit", None),
                ("plot_format", "pdf"),
                ("sheet", None),
                ("plate_id", ""),
            )
        },
        "fia": (
            {
                str(path): hashlib.sha256(data).hexdigest()
                for path, data in fia["files"].items()
            }
            if fia
            else None
        ),
    }


def verify_selection(selection, paths):
    actual = {key: sha256_file(path) for key, path in paths.items()}
    if actual != selection["identity"]["inputs"]:
        raise ValidationError(
            "A planned report input changed before processing"
        )
    state = report_inputs.read_collection_status(selection["path"])
    if (
        state.get("selected", {}) != selection["identity"]["selected_runs"]
        or state.get("processing_exclusions", [])
        != selection["identity"]["exclusions"]
    ):
        raise ValidationError(
            "Collection manifest changed after report planning"
        )
    if selection["fia"]:
        for path, data in selection["fia"]["files"].items():
            if path.read_bytes() != data:
                raise ValidationError(
                    f"FIA input changed after planning: {path}"
                )


def save_artifacts(directory):
    """Record outputs; logs and the completion marker remain mutable."""
    records = []
    for path in sorted(directory.rglob("*")):
        if not path.is_file() or path.is_symlink():
            continue
        if path.name in {
            "run_status.json",
            "run.log",
            "run_log.csv",
            "artifacts.json",
        }:
            continue
        records.append(
            {
                "path": str(path.relative_to(directory)),
                "sha256": sha256_file(path),
            }
        )
    save_json(directory / "artifacts.json", records)
    return sha256_file(directory / "artifacts.json")


def verified_previous(parent, fingerprint, log):
    """A folder alone is insufficient: verify the completed output files."""
    for folder in sorted(parent.glob("UMA_Report_*"), reverse=True):
        if folder.is_symlink() or not folder.is_dir():
            continue
        marker = folder / "run_status.json"
        try:
            if marker.is_symlink():
                continue
            state = json.loads(marker.read_text())
            if (
                state.get("status") not in {"SUCCESS", "PARTIAL"}
                or state.get("report_fingerprint") != fingerprint
            ):
                continue
            manifest = inside(folder, "artifacts.json")
            if sha256_file(manifest) != state.get("artifacts_sha256"):
                raise ValidationError("Output manifest changed")
            records = json.loads(manifest.read_text())
            if not records or not any(
                r["path"].endswith(".xlsx") for r in records
            ):
                raise ValidationError("No workbook in output manifest")
            for record in records:
                path = inside(folder, record["path"])
                if sha256_file(path) != record["sha256"]:
                    raise ValidationError(f"Output changed: {path}")
            return folder
        except (OSError, ValueError, KeyError, TypeError) as error:
            log.event("WARNING", "Previous report", f"{folder}: {error}")
    return None


def choose(plans, args):
    mode = getattr(args, "report_mode", None)
    requested = getattr(args, "area_variants", None)
    if requested:
        mode = "select"
    if mode is None and getattr(args, "wizard", False):
        mode = _ask(
            "Area reports: changed, all, or select",
            "changed",
            _choice({key: key for key in ("changed", "all", "select")}),
        )
    mode = mode or "changed"
    if mode == "select" and not requested:

        def selected_numbers(value):
            values = [int(x.strip()) for x in value.split(",")]
            if not values or any(x < 1 or x > len(plans) for x in values):
                raise ValueError(
                    "Choose variant numbers from the displayed list"
                )
            return [plans[x - 1]["variant"]["id"] for x in values]

        requested = _ask(
            "Variant numbers, separated by commas", "1", selected_numbers
        )
    if mode == "select":
        known = {item["variant"]["id"] for item in plans}
        unknown = set(requested) - known
        if unknown:
            raise ValidationError(
                "Unknown Area variant IDs: " + ", ".join(sorted(unknown))
            )
    return [
        item
        for item in plans
        if mode == "all"
        or (mode == "changed" and item.get("previous") is None)
        or (mode == "select" and item["variant"]["id"] in requested)
    ]


@folder_logged("source")
def process_source(source, input_json, args):
    from . import report

    buffered = report.BufferedLog()
    try:
        collection = report_inputs.select_combined(source, buffered)
        state = report_inputs.read_collection_status(collection)
        if state.get("collection_schema") != COLLECTION_SCHEMA:
            return report.process_folder.__wrapped__(source, input_json, args)
        variants = state.get("variants")
        if not isinstance(variants, list) or not variants:
            raise ValidationError("Collection has no Area variants")
        ids = [variant.get("id") for variant in variants]
        if len(set(ids)) != len(ids):
            raise ValidationError("Duplicate Area variant IDs in collection")
        collection_errors = state.get("errors", [])
        for error in collection_errors:
            buffered.event("WARNING", "Collection incomplete", str(error))
        fia = None
        if getattr(args, "with_fia_nuclei", False):
            from .report_nuclei import select_run

            fia = select_run(source, buffered)
        plans = []
        for index, variant in enumerate(variants, 1):
            item = {"variant": variant, "collection": collection, "fia": fia}
            plans.append(item)
            try:
                if variant.get("status") != "SUCCESS":
                    raise ValidationError(
                        f"Variant collection is {variant.get('status')}: "
                        + "; ".join(variant.get("errors", []))
                    )
                path = inside(collection, variant["path"])
                paths = report_inputs.discover_inputs(
                    path, source.name, collection
                )
                variant_state = report_inputs.read_collection_status(path)
                if sha256_file(inside(path, "processing_exclusions.csv")) != (
                    variant_state.get("exclusions_sha256")
                ):
                    raise ValidationError(
                        "Processing exclusions are missing or changed"
                    )
                with paths["fibronectin"].open(
                    encoding="utf-8-sig", newline=""
                ) as stream:
                    actual_profile = profile(csv.DictReader(stream))
                if (
                    actual_profile != variant["profile"]
                    or actual_profile != variant_state.get("area_profile")
                    or variant_id(actual_profile) != variant["id"]
                    or variant_label(actual_profile) != variant["label"]
                ):
                    raise ValidationError(
                        "Area profile in collection disagrees with its CSV"
                    )
                resolved = identity(
                    source, args, paths, variant_state, variant, fia
                )
                signature = digest(resolved)
                previous = verified_previous(
                    assay_directory(source), signature, buffered
                )
                item.update(
                    path=path,
                    identity=resolved,
                    fingerprint=signature,
                    previous=previous,
                    previous_status=(
                        json.loads((previous / "run_status.json").read_text())[
                            "status"
                        ]
                        if previous
                        else None
                    ),
                )
                action = "UNCHANGED" if previous else "NEW / CHANGED"
                console(
                    f"{index}. {variant['label']} [{variant['id']}]: {action}"
                )
            except (
                OSError,
                ValueError,
                KeyError,
                TypeError,
                ReportValidationError,
            ) as error:
                item["error"] = str(error)
                console(
                    f"{index}. {variant.get('label')}: UNAVAILABLE: {error}"
                )
        chosen = choose(plans, args)
        # Planning itself is auditable, even if every output can be reused.
        from .run import RunLog, unique_output

        _, batch = unique_output(assay_directory(source), "Report_Batch_")
        log = RunLog(batch)
        try:
            buffered.replay(log)
            chosen_ids = {item["variant"]["id"] for item in chosen}
            rows = [
                {
                    "Variant_ID": item["variant"]["id"],
                    "Area_variant": item["variant"]["label"],
                    "Action": "PROCESS"
                    if item["variant"]["id"] in chosen_ids
                    else "REUSE"
                    if item.get("previous")
                    else "NOT_SELECTED",
                    "Previous_report": str(item.get("previous") or ""),
                    "Reason": item.get(
                        "error",
                        "Inputs/settings unchanged"
                        if item.get("previous")
                        else "New, changed, or missing output",
                    ),
                }
                for item in plans
            ]
            save_csv(batch / "report_plan.csv", PLAN_COLUMNS, rows)
            save_json(
                batch / "run_status.json",
                {
                    "status": "RUNNING",
                    "collection": str(collection),
                    "source_folder": str(source),
                    "plan": rows,
                },
            )
            for row in rows:
                log.event(
                    "INFO",
                    "Report plan",
                    f"{row['Area_variant']}: {row['Action']}; {row['Reason']}",
                )
            results = []
            for index, item in enumerate(chosen, 1):
                phase(
                    f"Report {index}/{len(chosen)}: {item['variant']['label']}"
                )
                if "error" in item:
                    log.event("ERROR", "Variant unavailable", item["error"])
                    results.append(report.ReportStatus.VALIDATION_FAILED)
                    continue
                result = report.process_folder.__wrapped__(
                    source, input_json, args, item
                )
                results.append(result)
                log.event(
                    result.value, "Variant finished", item["variant"]["label"]
                )
            bad = any(r is not report.ReportStatus.SUCCESS for r in results)
            unavailable = any("error" in item for item in plans)
            reused_partial = any(
                item.get("previous_status") == "PARTIAL"
                and item["variant"]["id"] not in chosen_ids
                for item in plans
            )
            if bad or unavailable or reused_partial or collection_errors:
                final = report.ReportStatus.PARTIAL
            else:
                final = report.ReportStatus.SUCCESS
            save_json(
                batch / "run_status.json",
                {
                    "status": final.value,
                    "source_folder": str(source),
                    "collection": str(collection),
                    "plan": rows,
                    "results": [r.value for r in results],
                    "collection_errors": collection_errors,
                },
            )
            outcome(
                final.value,
                f"{source}: {len(chosen)} variant(s) processed; "
                f"{len(plans) - len(chosen)} not recalculated. "
                f"Batch log: {batch}",
            )
            return final
        except BaseException as error:
            save_json(
                batch / "run_status.json",
                {
                    "status": "CANCELLED"
                    if isinstance(error, (KeyboardInterrupt, EOFError))
                    else "ERROR",
                    "source_folder": str(source),
                    "collection": str(collection),
                    "error": str(error),
                },
            )
            raise
        finally:
            log.close()
    except (KeyboardInterrupt, EOFError):
        raise KeyboardInterrupt
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        ReportValidationError,
    ) as error:
        return report.diagnostic_failure(source, input_json, error, buffered)
