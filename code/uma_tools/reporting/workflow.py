#!/usr/bin/env python3
"""
Generate one validated UMA report per source folder in input_paths.json.

Run: uma_report -i input_paths.json --fn-threshold 20 The latest
completed Combined_Results is selected before checking its three
summaries and manually supplied plate map. An incomplete selected input
does not cause a silent fallback to another dataset. No ImageJ runtime
is started.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from ..common.config import read_config
from ..common.files import safe_label, save_json
from ..common.version import package_version
from . import engine
from .collection import (
    COMBINED_PATTERN as COMBINED_PATTERN,
)
from .collection import (
    ROLES as ROLES,
)
from .collection import (
    SUCCESS_STATES as SUCCESS_STATES,
)
from .collection import (
    archive_inputs as archive_inputs,
)
from .collection import (
    collection_entries as collection_entries,
)
from .collection import (
    discover_inputs as discover_inputs,
)
from .collection import (
    read_collection_status as read_collection_status,
)
from .collection import (
    recorded_basename as recorded_basename,
)
from .collection import (
    select_combined as select_combined,
)
from .collection import (
    verify_sources as verify_sources,
)

SCRIPT_VERSION = "4.0.0"
ValidationError = engine.ValidationError


class BufferedLog:
    """
    Hold read-only discovery events until the output location is known.
    """

    def __init__(self):
        self.events = []

    def event(self, level, stage, message, console=True, timestamp=None):
        self.events.append(
            (
                level,
                stage,
                str(message),
                console,
                timestamp or engine.utc_now(),
            )
        )

    def replay(self, log):
        for level, stage, message, console, timestamp in self.events:
            log.event(
                level, stage, message, console=console, timestamp=timestamp
            )


def new_run(parent, source_name):
    label = safe_label(source_name)
    stem = (
        datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
        + f"_{os.getpid()}"
    )
    for counter in range(10000):
        run_id = stem + (f"_{counter:03d}" if counter else "")
        directory = parent / f"UMA_Report_{label}_{run_id}"
        try:
            directory.mkdir()
            return run_id, directory
        except FileExistsError:
            continue
    raise OSError("Could not allocate a unique report directory")


def best_effort(description, action, *args, **kwargs):
    """
    Keep a secondary diagnostics failure from masking the original
    error.
    """
    try:
        action(*args, **kwargs)
    except Exception as error:
        try:
            print(
                f"Could not {description}: {error}",
                file=sys.stderr,
                flush=True,
            )
        except OSError:
            pass


def diagnostic_failure(source, input_json, error, buffered=None):
    """
    Use the source folder or CWD when no selected collection can hold
    logs.
    """
    parents = [source, Path.cwd()] if source is not None else [Path.cwd()]
    for parent in dict.fromkeys(parents):
        log = None
        try:
            run_id, directory = new_run(
                parent,
                source.name if source is not None else "configuration_error",
            )
            log = engine.RunLog(directory)
            if buffered is not None:
                buffered.replay(log)
            log.event("FAILED", "Startup", str(error))
            save_json(
                directory / "run_status.json",
                {
                    "run_id": run_id,
                    "script_version": SCRIPT_VERSION,
                    "package_version": package_version(),
                    "status": "VALIDATION_FAILED"
                    if isinstance(error, (engine.ValidationError, ValueError))
                    else "ERROR",
                    "ended_utc": engine.utc_now(),
                    "source_folder": str(source)
                    if source is not None
                    else None,
                    "input_json": str(input_json),
                    "error": str(error),
                    "run_directory": str(directory),
                },
            )
            engine.save_details(
                directory / "validation_errors.csv",
                [{"Stage": "Startup", "Issue": str(error)}],
            )
            log.event("INFO", "Diagnostics", directory)
            return False, directory
        except OSError:
            continue
        finally:
            if log is not None:
                best_effort("close startup logs", log.close)
    print(
        f"Could not start report or save diagnostics: {error}", file=sys.stderr
    )
    return False, None


def process_folder(source, input_json, args):
    buffered = BufferedLog()
    try:
        combined = select_combined(source, buffered)
        run_id, directory = new_run(combined, source.name)
    except (OSError, ValueError, engine.ValidationError) as error:
        return diagnostic_failure(source, input_json, error, buffered)
    log = engine.RunLog(directory)
    status = {
        "run_id": run_id,
        "script_version": SCRIPT_VERSION,
        "package_version": package_version(),
        "status": "RUNNING",
        "stage": "Initialization",
        "started_utc": engine.utc_now(),
        "source_folder": str(source),
        "source_name": source.name,
        "combined_results_folder": str(combined),
        "input_json": str(input_json),
        "run_directory": str(directory),
    }
    candidate, final_path = None, None

    def stage(name):
        status["stage"] = name
        save_json(directory / "run_status.json", status)

    try:
        stage("Initialization")
        log.event(
            "STARTED",
            "Run",
            f"UMA-tools {package_version()}, report {SCRIPT_VERSION}; "
            f"source {source}",
        )
        log.event("INFO", "Output", directory)
        buffered.replay(log)
        parameters = {
            **status,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "python_executable": sys.executable,
            "fn_area_threshold_percent": args.fn_threshold,
            "template_sheet": args.sheet or "First worksheet",
            "plate_id": args.plate_id or source.name,
            "selection_rule": (
                "Latest completed collection; no fallback for missing/invalid "
                "selected inputs"
            ),
            "image_id_rule": (
                "Complete original filename including extension; exact "
                "alignment suffix mapping"
            ),
            "fn_filter_rule": (
                "FN_Area_Percent < threshold is excluded from filtered views; "
                "equality is retained"
            ),
            "filter_scope": (
                "Individual images; all full-data rows remain available"
            ),
            "replicates": (
                "Image-level points; wells are technical replicates; "
                "biological replicate is blank"
            ),
            "statistical_tests": "Not performed",
            "quartiles": "Linear interpolation (R type 7)",
            "whiskers": "1.5 x IQR",
            "thickness_units": dict(engine.THICKNESS_UNITS),
            "alignment_axis": [0, 100],
            "fibronectin_axis": [0, 100],
            "thickness_axis": (
                "Starts at zero; full-data limits shared by full/filtered "
                "pairs"
            ),
            "empty_groups": (
                "Preserve X positions with n=0; n=1 displays one point "
                "without a box"
            ),
        }
        save_json(directory / "run_parameters.json", parameters)
        stage("Parameter validation")
        threshold = engine.validate_fn_threshold(args.fn_threshold)
        stage("Input discovery")
        paths = discover_inputs(combined, source.name)
        status["input_files"] = {
            role: str(path) for role, path in paths.items()
        }
        log.event(
            "PASS",
            "Inputs",
            "Three unchanged summaries and one plate-template workbook found",
        )
        stage("Input archive")
        snapshots, manifests = archive_inputs(
            paths, input_json, combined, directory
        )
        parameters.update(
            input_files=status["input_files"],
            archived_inputs={
                role: str(path) for role, path in snapshots.items()
            },
        )
        stage("Dependencies")
        parameters["dependency_versions"] = engine.load_dependencies(log)
        save_json(directory / "run_parameters.json", parameters)
        stage("Input validation")
        data = engine.validate_and_merge(
            snapshots,
            args.sheet,
            args.plate_id or source.name,
            directory,
            log,
            threshold,
        )
        parameters.update(
            alignment_metric=data["metric"],
            alignment_angle=data["angle_label"],
            template_sheet=data["template_sheet"],
            plate_id=data["plate_id"],
            group_order=data["group_order"],
            source_field_map=data["field_map"],
            full_data_images=len(data["rows"]),
            filtered_images=len(data["retained_rows"]),
            excluded_from_filtered_plots=len(data["excluded_rows"]),
        )
        save_json(directory / "run_parameters.json", parameters)
        for filename, rows in (
            ("merged_data.csv", data["rows"]),
            ("filtered_data.csv", data["retained_rows"]),
            ("excluded_data.csv", data["excluded_rows"]),
        ):
            engine.save_csv(directory / filename, data["columns"], rows)
        for filename, key in (
            ("group_filter_counts.csv", "group_filter_counts"),
            ("well_filter_counts.csv", "well_filter_counts"),
        ):
            engine.save_csv(
                directory / filename, list(data[key][0]), data[key]
            )
        engine.save_csv(
            directory / "qc.csv", ["Check", "Value", "Details"], data["qc"]
        )
        status.update(
            total_images=len(data["rows"]),
            full_data_images=len(data["rows"]),
            included_images=len(data["retained_rows"]),
            excluded_images=len(data["excluded_rows"]),
            groups=len(data["group_order"]),
            wells=len(data["well_counts"]),
            annotation_coverage_percent=100,
        )
        stage("Plots")
        plots = engine.create_plots(data, directory / "Plots", log)
        save_json(directory / "plot_manifest.json", plots)
        stage("Workbook export")
        candidate = directory / "report_pending.xlsx"
        final_path = (
            directory / f"UMA_Report_{safe_label(source.name)}_{run_id}.xlsx"
        )
        workbook = engine.build_workbook(data, plots, log.events, run_id)
        try:
            workbook.save(candidate)
            engine.verify_workbook(candidate, data)
            log.event(
                "PASS",
                "Workbook verification",
                "All data, 20 sheets, and 13 embedded plots verified",
            )
            completion_time = engine.utc_now()
            message = (
                f"{len(data['rows'])} images; "
                f"{len(data['retained_rows'])} retained; "
                f"{len(data['excluded_rows'])} excluded from filtered views; "
                f"13 plots. Workbook: {final_path.name}"
            )
            completion = dict(
                zip(
                    engine.EVENT_COLUMNS,
                    [completion_time, "SUCCESS", "Run", message],
                )
            )
            log_sheet = workbook["Run Log"]
            log_sheet.delete_rows(1, log_sheet.max_row)
            engine.write_table(
                log_sheet,
                engine.EVENT_COLUMNS,
                log.events + [completion],
                widths=[28, 14, 32, 130],
            )
            workbook.save(candidate)
        finally:
            workbook.close()
        engine.verify_workbook(candidate, data, require_success=True)
        verify_sources(manifests)
        candidate.rename(final_path)
        log.event("SUCCESS", "Run", message, timestamp=completion_time)
        status.update(
            status="SUCCESS",
            stage="Completed",
            ended_utc=completion_time,
            workbook=str(final_path),
            generated_plots=13,
        )
        save_json(directory / "run_status.json", status)
        print(f"Workbook: {final_path}\nLog: {log.path}", flush=True)
        return True, directory
    except (Exception, KeyboardInterrupt) as error:
        validation = isinstance(error, engine.ValidationError)
        label = (
            "VALIDATION_FAILED"
            if validation
            else "CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else "ERROR"
        )
        message = str(error) or "Run interrupted by the user"
        trace = traceback.format_exc()
        # Cleanup and failure status must not depend on successful log
        # writes.
        for path in (candidate, final_path):
            if path is not None:
                best_effort(
                    f"remove incomplete workbook {path}",
                    path.unlink,
                    missing_ok=True,
                )
        status.update(status=label, ended_utc=engine.utc_now(), error=message)
        status.pop("workbook", None)
        best_effort(
            "save failure status",
            save_json,
            directory / "run_status.json",
            status,
        )
        best_effort(
            "record failure in run log",
            log.event,
            "FAILED",
            status["stage"],
            message,
        )
        details = (
            error.details
            if validation and error.details
            else [{"Stage": status["stage"], "Issue": message}]
        )
        best_effort(
            "save diagnostic details",
            engine.save_details,
            directory
            / ("validation_errors.csv" if validation else "error_details.csv"),
            details,
        )
        for detail in details:
            best_effort(
                "record diagnostic detail",
                log.event,
                "ERROR",
                "Diagnostic detail",
                json.dumps(detail, ensure_ascii=False),
                console=False,
            )
        best_effort(
            "save traceback",
            (directory / "traceback.txt").write_text,
            trace,
            encoding="utf-8",
        )
        print(
            f"Report failed. Diagnostics: {directory}",
            file=sys.stderr,
            flush=True,
        )
        if isinstance(error, KeyboardInterrupt):
            raise
        return False, directory
    finally:
        best_effort("close report logs", log.close)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Create UMA plots and Excel reports from collected assay results"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="JSON file containing folder_paths",
    )
    parser.add_argument(
        "--fn-threshold",
        type=float,
        default=20.0,
        help=(
            "FN coverage cutoff in percent, 0-100 (default: 20); equality is "
            "retained"
        ),
    )
    parser.add_argument(
        "--sheet",
        default=None,
        help="Exact plate-template worksheet name (default: first sheet)",
    )
    parser.add_argument(
        "--plate-id",
        default="",
        help="Plate label (default: each source folder's name)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {package_version()} (report {SCRIPT_VERSION})",
    )
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
    print(
        f"Reports finished: {succeeded} folder(s) succeeded; {failed} failed.",
        flush=True,
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
