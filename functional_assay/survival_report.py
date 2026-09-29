"""Command for one plate's multi-day distributions and baseline changes."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import traceback
from collections.abc import Sequence
from importlib.metadata import version
from pathlib import Path

from uma_tools.files import safe_label, save_csv, save_json, sha256_file
from uma_tools.report_schema import ValidationError
from uma_tools.run import RunLog, unique_output, utc_now

from . import report_data, survival_data
from .functional_report import verify_inputs


def parse_args(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(
        description=(
            "Report one plate over explicit days and "
            "compare well changes with control"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Survival-assay JSON configuration",
    )
    parser.add_argument(
        "--stats-unit",
        choices=("well",),
        default=None,
        help="Enable Welch + Holm on well changes; omit for descriptive plots",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('uma-functional-assay')}",
    )
    return parser.parse_args(argv)


def archive_input(source, relative, output, manifest, role):
    """Keep separate snapshots for equal filenames from different days."""
    if source.is_symlink() or not source.is_file():
        raise ValidationError(f"Missing regular {role} input: {source}")
    digest = sha256_file(source)
    destination = output / "inputs" / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    if sha256_file(destination) != digest or sha256_file(source) != digest:
        raise ValidationError(f"Input changed while being archived: {source}")
    manifest.append(
        {
            "Role": role,
            "Source": str(source),
            "Snapshot": str(destination.relative_to(output)),
            "SHA256": digest,
        }
    )
    save_json(output / "input_manifest.json", manifest, allow_nan=False)
    return destination


def load_days(config, output, log, manifest):
    measurements, selections = {}, []
    for point in config["timepoints"]:
        day, folder = point["day"], point["folder"]
        messages = []
        if not folder.is_dir():
            raise ValidationError(f"Day {day}: folder not found: {folder}")
        try:
            selected = report_data.select_analysis(folder, messages)
        finally:
            for message in messages:
                log.event(
                    "WARNING", "Analysis selection", f"Day {day}: {message}"
                )
        log.event("INFO", "Analysis selection", f"Day {day}: {selected}")
        snapshot_dir = Path(f"day_{day}")
        status_path = archive_input(
            selected / "run_status.json",
            snapshot_dir / "run_status.json",
            output,
            manifest,
            f"Day {day} completion record",
        )
        summary_path = archive_input(
            selected / "Cell_Analysis_Summary.csv",
            snapshot_dir / "Cell_Analysis_Summary.csv",
            output,
            manifest,
            f"Day {day} measurements",
        )
        status = report_data.completed_status(status_path)
        measurements[day] = report_data.read_measurements(summary_path, status)
        selections.append(
            {
                "Day": day,
                "Source_Folder": str(folder),
                "Selected_Analysis": str(selected),
                "Summary_File": str(selected / "Cell_Analysis_Summary.csv"),
                "Completed_Wells": len(measurements[day]),
                "Parameters": json.dumps(status["parameters"], sort_keys=True),
            }
        )
    return measurements, selections


def run_report(config, input_path: Path, statistics: bool, input_digest=None):
    """Write a new experiment report without pooling days as replicates."""
    parent = config["output_dir"]
    parent.mkdir(parents=True, exist_ok=True)
    run_id, output = unique_output(
        parent, f"Survival_Report_{safe_label(config['experiment_name'])}_"
    )
    log = RunLog(output)
    status_path = output / "run_status.json"
    status = {
        "status": "running",
        "started_utc": utc_now(),
        "experiment_name": config["experiment_name"],
        "output": str(output),
        "functional_assay_version": version("uma-functional-assay"),
        "uma_tools_version": version("uma-tools"),
        "baseline_day": config["baseline_day"],
        "difference_days": config["difference_days"],
        "stats_unit": "well" if statistics else None,
    }
    pending = output / "report_pending.xlsx"
    manifest = []
    try:
        save_json(status_path, status, allow_nan=False)
        log.event("STARTED", "Survival report", config["experiment_name"])
        log.event("INFO", "Output", output)
        archive_input(
            input_path,
            Path("survival_config.json"),
            output,
            manifest,
            "Survival configuration",
        )
        if input_digest is not None and manifest[0]["SHA256"] != input_digest:
            raise ValidationError(
                "Survival configuration changed after it was read"
            )
        template = archive_input(
            config["plate_template"],
            Path("plate_template.xlsx"),
            output,
            manifest,
            "Plate template",
        )
        # Unmapped wells are warnings here: do not pass measured_wells.
        plate = report_data.read_design(
            template, config["sheet"], statistics=statistics
        )
        measurements, selections = load_days(config, output, log, manifest)
        rows, coverage, warnings = survival_data.join_days(
            measurements, plate, config["baseline_day"]
        )
        if not any(row["Annotation_Status"] == "MAPPED" for row in rows):
            raise ValidationError(
                "No measured wells match the plate-map annotations"
            )
        warnings = plate["warnings"] + warnings
        for warning in warnings:
            log.event("WARNING", "Well/parameter checks", warning)
        changes = survival_data.calculate_changes(
            rows, plate, config["baseline_day"], config["difference_days"]
        )
        for day in config["difference_days"]:
            eligible = {
                r["Well"]
                for r in changes
                if r["Day"] == day and r["Status"] == "PAIRED"
            }
            log.event(
                "INFO",
                "Paired changes",
                f"{day}-{config['baseline_day']}: "
                f"{len(eligible)} matched wells",
            )
        data = {
            **plate,
            "rows": rows,
            "changes": changes,
            "coverage": coverage,
            "selections": selections,
            "warnings": warnings,
            "baseline_day": config["baseline_day"],
            "difference_days": config["difference_days"],
            "days": sorted(measurements),
            "statistics_enabled": statistics,
        }
        data["summary"] = survival_data.summarize(data)
        data["comparisons"] = (
            survival_data.compare_changes(data) if statistics else []
        )
        for comparison in data["comparisons"]:
            if comparison["Status"] == "Not tested":
                log.event(
                    "WARNING",
                    "Statistics",
                    f"{comparison['Comparison_Block']} "
                    f"{comparison['Comparison']}: {comparison['Treatment']} "
                    f"versus {comparison['Control']}, "
                    f"{comparison['Metric']}: {comparison['Reason']}",
                )
        log.event(
            "INFO",
            "Statistics",
            "Welch on changes; Holm across both endpoints "
            "and all requested days per color"
            if statistics
            else "Disabled; --stats-unit omitted",
        )
        from .survival_output import (
            build_workbook,
            render_plots,
            tables,
            verify_workbook,
        )

        for _, name, columns, records in tables(data):
            save_csv(
                output / name,
                columns,
                ({key: row.get(key) for key in columns} for row in records),
                encoding="utf-8-sig",
            )
        plots = render_plots(data, output, config["experiment_name"])
        build_workbook(data, plots, template, pending, status)
        verify_workbook(pending, data, plots)
        verify_inputs(manifest)
        name = safe_label(config["experiment_name"])
        workbook = output / f"Survival_Report_{name}_{run_id}.xlsx"
        pending.replace(workbook)
        save_json(
            output / "plot_manifest.json",
            [
                {
                    "View": plot["view"],
                    "Metric": plot["metric"],
                    "File": plot["path"].name,
                    "Points": plot["points"],
                    "Annotations": plot["annotations"],
                }
                for plot in plots
            ],
            allow_nan=False,
        )
        status.update(
            status="completed",
            workbook=workbook.name,
            days=data["days"],
            raw_measurements=len(rows),
            unmapped_measurements=sum(
                r["Annotation_Status"] == "UNMAPPED" for r in rows
            ),
            paired_well_changes=sum(r["Status"] == "PAIRED" for r in changes)
            // len(report_data.METRICS),
            planned_comparisons=len(data["comparisons"]),
            tested_comparisons=sum(
                r["Status"] == "Tested" for r in data["comparisons"]
            ),
            selections=selections,
            warnings=warnings,
            inputs=manifest,
            plots=[plot["path"].name for plot in plots],
        )
        log.event("SUCCESS", "Survival report", workbook)
    except (Exception, KeyboardInterrupt) as error:
        status.update(
            status="failed", error=str(error) or type(error).__name__
        )
        pending.unlink(missing_ok=True)
        log.event("FAILED", "Survival report", status["error"])
        log.event("ERROR", "Traceback", traceback.format_exc(), console=False)
        if isinstance(error, ValidationError) and error.details:
            columns = list(
                dict.fromkeys(key for row in error.details for key in row)
            )
            save_csv(
                output / "validation_errors.csv",
                columns,
                error.details,
                encoding="utf-8-sig",
            )
        if isinstance(error, KeyboardInterrupt):
            raise
    finally:
        status["finished_utc"] = utc_now()
        try:
            save_json(status_path, status, allow_nan=False)
        finally:
            log.close()
    return status


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    path = Path(args.input).expanduser().absolute()
    try:
        if path.name.startswith("._"):
            raise ValidationError("AppleDouble ._ JSON files are not inputs")
        digest = sha256_file(path)
        config = survival_data.read_config(path)
        if sha256_file(path) != digest:
            raise ValidationError(
                "Survival configuration changed while being read"
            )
        result = run_report(config, path, args.stats_unit is not None, digest)
    except KeyboardInterrupt:
        print("Survival report interrupted.", file=sys.stderr)
        return 130
    except Exception as error:
        print(f"Survival report failed: {error}", file=sys.stderr)
        return 1
    if result["status"] != "completed":
        print(
            f"Survival report failed. Diagnostics: {result['output']}",
            file=sys.stderr,
        )
        return 1
    print(f"Survival report completed: {result['output']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
