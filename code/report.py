#!/usr/bin/env python3
"""Generate one validated UMA report per source folder in input_paths.json.

Run: uma_report -i input_paths.json --fn-threshold 20
The latest completed Combined_Results is selected before checking its three
summaries and manually supplied plate map. An incomplete selected input does
not cause a silent fallback to another dataset. No ImageJ runtime is started.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path, PureWindowsPath
import platform
import re
import sys
import traceback

if __package__:
    from . import report_rendering as engine
    from .collect_results import read_config, safe_label, save_json
else:
    import report_rendering as engine
    from collect_results import read_config, safe_label, save_json


SCRIPT_VERSION = "4.0.0"
ValidationError = engine.ValidationError
ROLES = {"Alignment": ("alignment", "Alignment_Summary.csv"),
         "Thickness": ("thickness", "Thickness_Summary.csv"),
         "Area": ("fibronectin", "Fibronectin_Area_Summary.csv")}
COMBINED_PATTERN = re.compile(r"Combined_Results_.+_(\d{8}_\d{6})(?:_(\d{6}))?(?:_\d+)?")
SUCCESS_STATES = {"SUCCESS", "SUCCESS_WITH_MISSING_ANALYSES"}


def package_version():
    try:
        return version("uma-tools")
    except PackageNotFoundError:
        return "not installed"


class BufferedLog:
    """Hold read-only discovery events until the output location is known."""

    def __init__(self):
        self.events = []

    def event(self, level, stage, message, console=True, timestamp=None):
        self.events.append((level, stage, str(message), console, timestamp or engine.utc_now()))

    def replay(self, log):
        for level, stage, message, console, timestamp in self.events:
            log.event(level, stage, message, console=console, timestamp=timestamp)


def read_collection_status(combined):
    path = combined / "run_status.json"
    if path.is_symlink() or not path.is_file():
        raise engine.ValidationError(f"Collector completion record is missing or not a regular file: {path}")
    status = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(status, dict):
        raise engine.ValidationError(f"Invalid collector completion record: {path}")
    return status


def select_combined(source, log):
    """Choose the newest completed collection, without inspecting its template."""
    candidates = []
    for path in sorted(source.iterdir()):
        if path.name.startswith(".") or not path.name.startswith("Combined_Results_"):
            continue
        if path.is_symlink() or not path.is_dir():
            continue
        try:
            match = COMBINED_PATTERN.fullmatch(path.name)
            if not match:
                raise engine.ValidationError("Unrecognized Combined_Results timestamp")
            timestamp = datetime.strptime(match[1], "%Y%m%d_%H%M%S")
            if match[2]:
                timestamp = timestamp.replace(microsecond=int(match[2]))
            status = read_collection_status(path)
            if status.get("status") not in SUCCESS_STATES:
                raise engine.ValidationError(f"Collector status is {status.get('status')!r}")
            candidates.append((timestamp, path))
        except (OSError, ValueError, engine.ValidationError) as error:
            log.event("WARNING", "Collection selection", f"Skipping {path}: {error}")
    if not candidates:
        raise engine.ValidationError(f"No successful Combined_Results directory found in {source}")
    candidates.sort(key=lambda item: (item[0], item[1].name), reverse=True)
    latest = [path for timestamp, path in candidates if timestamp == candidates[0][0]]
    if len(latest) != 1:
        raise engine.ValidationError("Several successful collections have the same latest timestamp: "
                                     + ", ".join(str(path) for path in latest))
    log.event("INFO", "Collection selection", f"Selected latest successful collection: {latest[0]}")
    return latest[0]


def recorded_basename(value):
    """Use archived filenames even after a collection is moved between drives."""
    if not isinstance(value, str) or not value:
        raise engine.ValidationError("Collector manifest contains an empty CSV path")
    name = PureWindowsPath(value).name if "\\" in value else Path(value).name
    if not name or name.startswith((".", "~$")) or name in {".", ".."}:
        raise engine.ValidationError(f"Invalid collected CSV filename: {value!r}")
    return name


def collection_entries(status):
    """Validate the complete collector manifest, including archived copies."""
    if not isinstance(status, dict):
        raise engine.ValidationError("Invalid collector completion record")
    if status.get("status") not in SUCCESS_STATES:
        raise engine.ValidationError("The selected collection is no longer marked successful")
    entries = status.get("copied_csvs")
    if not isinstance(entries, list):
        raise engine.ValidationError("Collector completion record has no copied_csvs manifest")
    by_analysis = {}
    for entry in entries:
        if not isinstance(entry, dict) or entry.get("analysis") not in ROLES:
            raise engine.ValidationError("Collector manifest contains an unknown analysis entry")
        analysis = entry["analysis"]
        if analysis in by_analysis:
            raise engine.ValidationError(f"Duplicate {analysis} entries in collector manifest")
        by_analysis[analysis] = entry
    missing = [analysis for analysis in ROLES if analysis not in by_analysis]
    if missing:
        raise engine.ValidationError("All three analyses are required. Missing: " + ", ".join(missing))
    return by_analysis


def discover_inputs(combined, source_name):
    """Require all three collected summaries and exactly one visible plate map."""
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
            issue = f"Collected {analysis} CSV is missing or not a regular file: {path}"
        elif not isinstance(digest, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
            issue = f"Missing or invalid {analysis} SHA256 in collector manifest"
        elif engine.sha256_file(path) != digest.lower():
            issue = f"Collected {analysis} CSV has changed since collection: {path}"
        if issue:
            details.append({"Input": role, "Path": str(path), "Issue": issue})
        else:
            paths[role] = path
    templates = sorted(path for path in combined.iterdir()
                       if path.suffix.lower() == ".xlsx" and not path.name.startswith((".", "~$"))
                       and path.is_file() and not path.is_symlink())
    if len(templates) != 1:
        issue = ("No plate-template .xlsx found in the selected Combined_Results directory"
                 if not templates else "More than one plate-template .xlsx found: " + ", ".join(path.name for path in templates))
        details.append({"Input": "template", "Path": str(combined), "Issue": issue})
    else:
        paths["template"] = templates[0]
    if details:
        raise engine.ValidationError("Selected collection is not ready for reporting. "
                                     + "; ".join(item["Issue"] for item in details), details)
    return paths


def new_run(parent, source_name):
    label = safe_label(source_name)
    stem = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f") + f"_{os.getpid()}"
    for counter in range(10000):
        run_id = stem + (f"_{counter:03d}" if counter else "")
        directory = parent / f"UMA_Report_{label}_{run_id}"
        try:
            directory.mkdir()
            return run_id, directory
        except FileExistsError:
            continue
    raise OSError("Could not allocate a unique report directory")


def archive_inputs(paths, input_json, combined, directory):
    """Make verified input snapshots before parsing or plotting their contents."""
    archive = directory / "Inputs"
    archive.mkdir()
    sources = {**paths, "configuration": input_json,
               "collector_status": combined / "run_status.json"}
    for role, name in (("collector_image_check", "image_check.csv"),
                       ("collector_selection_report", "selection_report.csv")):
        path = combined / name
        if path.is_file() and not path.is_symlink():
            sources[role] = path
    names = {"configuration": "input_paths.json", "collector_status": "collector_run_status.json",
             "collector_image_check": "collector_image_check.csv",
             "collector_selection_report": "collector_selection_report.csv"}
    records, snapshots, used = [], {}, set()
    for role, source in sources.items():
        name = names.get(role, source.name)
        if name in used:
            raise engine.ValidationError(f"Input archive filename collision: {name}")
        used.add(name)
        data = source.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        target = archive / name
        target.write_bytes(data)
        if engine.sha256_file(target) != digest:
            raise OSError(f"Input copy verification failed: {target}")
        records.append({"Input": role, "Path": str(source), "Archived_Path": str(target),
                        "Archive_Relative_Path": str(target.relative_to(directory)),
                        "Bytes": len(data), "SHA256": digest})
        if role in paths:
            snapshots[role] = target
    engine.save_csv(directory / "input_manifest.csv", list(records[0]), records)
    # Reconcile copied summary bytes against the archived collector record as
    # well, closing the gap between discovery and archiving.
    archived_status = json.loads((archive / "collector_run_status.json").read_text(encoding="utf-8-sig"))
    entries = collection_entries(archived_status)
    archived_records = {record["Input"]: record for record in records}
    for analysis, (role, _) in ROLES.items():
        entry, record = entries[analysis], archived_records[role]
        digest = entry.get("sha256")
        if (recorded_basename(entry.get("path")) != paths[role].name
                or not isinstance(digest, str) or digest.lower() != record["SHA256"]):
            raise engine.ValidationError("Collector manifest changed while inputs were copied: " + record["Path"])
    verify_sources(records)
    return snapshots, records


def verify_sources(records):
    for record in records:
        if engine.sha256_file(Path(record["Path"])) != record["SHA256"]:
            raise engine.ValidationError("An input changed during report generation: " + record["Path"])


def best_effort(description, action, *args, **kwargs):
    """Keep a secondary diagnostics failure from masking the original error."""
    try:
        action(*args, **kwargs)
    except Exception as error:
        try:
            print(f"Could not {description}: {error}", file=sys.stderr, flush=True)
        except OSError:
            pass


def diagnostic_failure(source, input_json, error, buffered=None):
    """Use the source folder or CWD when no selected collection can hold logs."""
    parents = [source, Path.cwd()] if source is not None else [Path.cwd()]
    for parent in dict.fromkeys(parents):
        log = None
        try:
            run_id, directory = new_run(parent, source.name if source is not None else "configuration_error")
            log = engine.RunLog(directory)
            if buffered is not None:
                buffered.replay(log)
            log.event("FAILED", "Startup", str(error))
            save_json(directory / "run_status.json", {
                "run_id": run_id, "script_version": SCRIPT_VERSION, "package_version": package_version(),
                "status": "VALIDATION_FAILED" if isinstance(error, (engine.ValidationError, ValueError)) else "ERROR",
                "ended_utc": engine.utc_now(), "source_folder": str(source) if source is not None else None,
                "input_json": str(input_json), "error": str(error), "run_directory": str(directory)})
            engine.save_details(directory / "validation_errors.csv", [{"Stage": "Startup", "Issue": str(error)}])
            log.event("INFO", "Diagnostics", directory)
            return False, directory
        except OSError:
            continue
        finally:
            if log is not None:
                best_effort("close startup logs", log.close)
    print(f"Could not start report or save diagnostics: {error}", file=sys.stderr)
    return False, None


def process_folder(source, input_json, args):
    buffered = BufferedLog()
    try:
        combined = select_combined(source, buffered)
        run_id, directory = new_run(combined, source.name)
    except (OSError, ValueError, engine.ValidationError) as error:
        return diagnostic_failure(source, input_json, error, buffered)
    log = engine.RunLog(directory)
    status = {"run_id": run_id, "script_version": SCRIPT_VERSION, "package_version": package_version(),
              "status": "RUNNING", "stage": "Initialization", "started_utc": engine.utc_now(),
              "source_folder": str(source), "source_name": source.name,
              "combined_results_folder": str(combined), "input_json": str(input_json),
              "run_directory": str(directory)}
    candidate, final_path = None, None

    def stage(name):
        status["stage"] = name
        save_json(directory / "run_status.json", status)

    try:
        stage("Initialization")
        log.event("STARTED", "Run", f"UMA-tools {package_version()}, report {SCRIPT_VERSION}; source {source}")
        log.event("INFO", "Output", directory)
        buffered.replay(log)
        parameters = {**status, "python": platform.python_version(), "platform": platform.platform(),
                      "python_executable": sys.executable, "fn_area_threshold_percent": args.fn_threshold,
                      "template_sheet": args.sheet or "First worksheet", "plate_id": args.plate_id or source.name,
                      "selection_rule": "Latest completed collection; no fallback for missing/invalid selected inputs",
                      "image_id_rule": "Complete original filename including extension; exact alignment suffix mapping",
                      "fn_filter_rule": "FN_Area_Percent < threshold is excluded from filtered views; equality is retained",
                      "filter_scope": "Individual images; all full-data rows remain available",
                      "replicates": "Image-level points; wells are technical replicates; biological replicate is blank",
                      "statistical_tests": "Not performed", "quartiles": "Linear interpolation (R type 7)",
                      "whiskers": "1.5 x IQR", "thickness_units": dict(engine.THICKNESS_UNITS),
                      "alignment_axis": [0, 100], "fibronectin_axis": [0, 100],
                      "thickness_axis": "Starts at zero; full-data limits shared by full/filtered pairs",
                      "empty_groups": "Preserve X positions with n=0; n=1 displays one point without a box"}
        save_json(directory / "run_parameters.json", parameters)
        stage("Parameter validation")
        threshold = engine.validate_fn_threshold(args.fn_threshold)
        stage("Input discovery")
        paths = discover_inputs(combined, source.name)
        status["input_files"] = {role: str(path) for role, path in paths.items()}
        log.event("PASS", "Inputs", "Three unchanged summaries and one plate-template workbook found")
        stage("Input archive")
        snapshots, manifests = archive_inputs(paths, input_json, combined, directory)
        parameters.update(input_files=status["input_files"],
                          archived_inputs={role: str(path) for role, path in snapshots.items()})
        stage("Dependencies")
        parameters["dependency_versions"] = engine.load_dependencies(log)
        save_json(directory / "run_parameters.json", parameters)
        stage("Input validation")
        data = engine.validate_and_merge(snapshots, args.sheet, args.plate_id or source.name,
                                         directory, log, threshold)
        parameters.update(alignment_metric=data["metric"], alignment_angle=data["angle_label"],
                          template_sheet=data["template_sheet"], plate_id=data["plate_id"],
                          group_order=data["group_order"], source_field_map=data["field_map"],
                          full_data_images=len(data["rows"]), filtered_images=len(data["retained_rows"]),
                          excluded_from_filtered_plots=len(data["excluded_rows"]))
        save_json(directory / "run_parameters.json", parameters)
        for filename, rows in (("merged_data.csv", data["rows"]), ("filtered_data.csv", data["retained_rows"]),
                               ("excluded_data.csv", data["excluded_rows"])):
            engine.save_csv(directory / filename, data["columns"], rows)
        for filename, key in (("group_filter_counts.csv", "group_filter_counts"),
                              ("well_filter_counts.csv", "well_filter_counts")):
            engine.save_csv(directory / filename, list(data[key][0]), data[key])
        engine.save_csv(directory / "qc.csv", ["Check", "Value", "Details"], data["qc"])
        status.update(total_images=len(data["rows"]), full_data_images=len(data["rows"]),
                      included_images=len(data["retained_rows"]), excluded_images=len(data["excluded_rows"]),
                      groups=len(data["group_order"]), wells=len(data["well_counts"]), annotation_coverage_percent=100)
        stage("Plots")
        plots = engine.create_plots(data, directory / "Plots", log)
        save_json(directory / "plot_manifest.json", plots)
        stage("Workbook export")
        candidate = directory / "report_pending.xlsx"
        final_path = directory / f"UMA_Report_{safe_label(source.name)}_{run_id}.xlsx"
        workbook = engine.build_workbook(data, plots, log.events, run_id)
        try:
            workbook.save(candidate)
            engine.verify_workbook(candidate, data)
            log.event("PASS", "Workbook verification", "All data, 20 sheets, and 13 embedded plots verified")
            completion_time = engine.utc_now()
            message = (f"{len(data['rows'])} images; {len(data['retained_rows'])} retained; "
                       f"{len(data['excluded_rows'])} excluded from filtered views; 13 plots. Workbook: {final_path.name}")
            completion = dict(zip(engine.EVENT_COLUMNS, [completion_time, "SUCCESS", "Run", message]))
            log_sheet = workbook["Run Log"]
            log_sheet.delete_rows(1, log_sheet.max_row)
            engine.write_table(log_sheet, engine.EVENT_COLUMNS, log.events + [completion], widths=[28, 14, 32, 130])
            workbook.save(candidate)
        finally:
            workbook.close()
        engine.verify_workbook(candidate, data, require_success=True)
        verify_sources(manifests)
        candidate.rename(final_path)
        log.event("SUCCESS", "Run", message, timestamp=completion_time)
        status.update(status="SUCCESS", stage="Completed", ended_utc=completion_time,
                      workbook=str(final_path), generated_plots=13)
        save_json(directory / "run_status.json", status)
        print(f"Workbook: {final_path}\nLog: {log.path}", flush=True)
        return True, directory
    except (Exception, KeyboardInterrupt) as error:
        validation = isinstance(error, engine.ValidationError)
        label = "VALIDATION_FAILED" if validation else "CANCELLED" if isinstance(error, KeyboardInterrupt) else "ERROR"
        message = str(error) or "Run interrupted by the user"
        trace = traceback.format_exc()
        # Cleanup and failure status must not depend on successful log writes.
        for path in (candidate, final_path):
            if path is not None:
                best_effort(f"remove incomplete workbook {path}", path.unlink, missing_ok=True)
        status.update(status=label, ended_utc=engine.utc_now(), error=message)
        status.pop("workbook", None)
        best_effort("save failure status", save_json, directory / "run_status.json", status)
        best_effort("record failure in run log", log.event, "FAILED", status["stage"], message)
        details = error.details if validation and error.details else [{"Stage": status["stage"], "Issue": message}]
        best_effort("save diagnostic details", engine.save_details,
                    directory / ("validation_errors.csv" if validation else "error_details.csv"), details)
        for detail in details:
            best_effort("record diagnostic detail", log.event, "ERROR", "Diagnostic detail",
                        json.dumps(detail, ensure_ascii=False), console=False)
        best_effort("save traceback", (directory / "traceback.txt").write_text, trace, encoding="utf-8")
        print(f"Report failed. Diagnostics: {directory}", file=sys.stderr, flush=True)
        if isinstance(error, KeyboardInterrupt):
            raise
        return False, directory
    finally:
        best_effort("close report logs", log.close)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Create UMA plots and Excel reports from collected assay results")
    parser.add_argument("-i", "--input", required=True, help="JSON file containing folder_paths")
    parser.add_argument("--fn-threshold", type=float, default=20.0,
                        help="FN coverage cutoff in percent, 0-100 (default: 20); equality is retained")
    parser.add_argument("--sheet", default=None, help="Exact plate-template worksheet name (default: first sheet)")
    parser.add_argument("--plate-id", default="", help="Plate label (default: each source folder's name)")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {package_version()} (report {SCRIPT_VERSION})")
    args = parser.parse_args(argv)
    input_json = Path(args.input).expanduser().absolute()
    try:
        folders = read_config(input_json)
    except (OSError, ValueError) as error:
        diagnostic_failure(None, input_json, error)
        return 1
    seen, succeeded, failed = set(), 0, 0
    for source in folders:
        try:
            canonical = source.resolve()
            if canonical in seen:
                print(f"Skipping duplicate JSON folder: {source}", flush=True)
                continue
            seen.add(canonical)
            success, _ = process_folder(source, input_json, args)
        except KeyboardInterrupt:
            print("Report generation interrupted.", file=sys.stderr)
            return 130
        except Exception as error:
            diagnostic_failure(source, input_json, error)
            success = False
        succeeded += int(success)
        failed += int(not success)
    print(f"Reports finished: {succeeded} folder(s) succeeded; {failed} failed.", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
