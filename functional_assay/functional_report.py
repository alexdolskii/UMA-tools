"""Create per-plate functional-assay reports from completed cell analysis."""

from __future__ import annotations

import argparse
import json
import shutil
import traceback
from collections.abc import Sequence
from importlib.metadata import version
from pathlib import Path

from uma_tools.config import read_config
from uma_tools.files import safe_label, save_csv, save_json, sha256_file
from uma_tools.plate_order import GROUP_ORDER_COLUMNS
from uma_tools.report_schema import ValidationError
from uma_tools.run import unique_output, utc_now

from . import report_data
from .plot_palette import palette_tables, save_palette
from .workflow import (
    NoInputError,
    RunLog,
    assay_directory,
    command_error,
    exclusions,
    save_exclusions,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot object count and mask area from the latest "
            "finalized cell analysis (including partial results)"
        )
    )
    parser.add_argument(
        "-i", "--input", required=True, help="UMA folder_paths JSON"
    )
    parser.add_argument(
        "--stats-unit",
        choices=("well",),
        default=None,
        help=(
            "Enable Welch + Holm control comparisons using wells; "
            "omit to disable statistics"
        ),
    )
    parser.add_argument(
        "--sheet", help="Plate-map worksheet name (default: first sheet)"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('uma-functional-assay')}",
    )
    return parser.parse_args(argv)


def snapshot_inputs(paths: dict, output: Path) -> tuple[dict, list]:
    """Archive exact inputs and detect writes occurring during the copy."""
    directory = output / "inputs"
    directory.mkdir()
    names = {
        "summary": "Cell_Analysis_Summary.csv",
        "analysis_status": "analysis_run_status.json",
        "input_json": "input_config.json",
    }
    snapshots, manifest = {}, []
    for role, source in paths.items():
        digest = sha256_file(source)
        destination = directory / names.get(role, source.name)
        shutil.copy2(source, destination)
        if sha256_file(destination) != digest or sha256_file(source) != digest:
            raise ValidationError(
                f"Input changed while being archived: {source}"
            )
        snapshots[role] = destination
        manifest.append(
            {
                "Role": role,
                "Source": str(source),
                "Snapshot": str(destination.relative_to(output)),
                "SHA256": digest,
            }
        )
    save_json(output / "input_manifest.json", manifest, allow_nan=False)
    return snapshots, manifest


def verify_inputs(manifest: list) -> None:
    for record in manifest:
        if sha256_file(Path(record["Source"])) != record["SHA256"]:
            raise ValidationError(
                f"Input changed during report generation: {record['Source']}"
            )


def _save_tables(data: dict, output: Path):
    from .report_output import WELL_COLUMNS

    tables = [
        (
            "Group_Order.csv",
            GROUP_ORDER_COLUMNS,
            data.get("group_order_records", []),
        ),
        ("Well_Data.csv", WELL_COLUMNS, data["rows"]),
        ("Condition_Summary.csv", report_data.GROUP_COLUMNS, data["summary"]),
        (
            "Plate_Coverage.csv",
            report_data.ANNOTATION_COLUMNS,
            data["diagnostics"],
        ),
    ]
    if data["statistics_enabled"]:
        tables.append(
            (
                "Statistics.csv",
                report_data.COMPARISON_COLUMNS,
                data["comparisons"],
            )
        )
    tables.extend(
        (name.replace(" ", "_") + ".csv", columns, rows)
        for name, columns, rows in palette_tables(data)
    )
    for name, columns, rows in tables:
        save_csv(
            output / name,
            columns,
            ({key: row.get(key) for key in columns} for row in rows),
            encoding="utf-8-sig",
        )


def process_folder(source: Path, input_json: Path, args) -> dict:
    """Keep one run per plate/time point, including diagnostics on failure."""
    selected, selection_error = None, None
    messages = []
    try:
        selected = report_data.select_analysis(source, messages)
    except Exception as error:
        selection_error = error
    run_id, output = unique_output(
        assay_directory(source),
        f"Functional_Report_{safe_label(source.name)}_",
    )
    log = RunLog(output, source, "functional_report")
    status_path = output / "run_status.json"
    status = {
        "status": "RUNNING",
        "started_utc": utc_now(),
        "source": str(source),
        "analysis": str(selected) if selected else None,
        "output": str(output),
        "functional_assay_version": version("uma-functional-assay"),
        "uma_tools_version": version("uma-tools"),
        "stats_unit": args.stats_unit,
    }
    pending = output / "report_pending.xlsx"
    try:
        save_json(status_path, status, allow_nan=False)
        log.event(
            "STARTED",
            "Report",
            f"Functional assay {status['functional_assay_version']}; {source}",
        )
        log.event("INFO", "Output", output)
        log.event(
            "INFO",
            "Plan",
            "1/4 inputs; 2/4 tables/statistics; "
            "3/4 two plots; 4/4 verify workbook",
        )
        log.phase(1, 4, "Validate inputs")
        for message in messages:
            log.event("WARNING", "Analysis selection", message)
        if selection_error is not None:
            raise selection_error
        log.event(
            "INFO",
            "Analysis selection",
            f"Latest finalized run: {selected}",
        )
        paths = report_data.discover_inputs(selected)
        paths["input_json"] = input_json
        snapshots, manifest = snapshot_inputs(paths, output)
        original_status = report_data.completed_status(
            snapshots["analysis_status"]
        )
        rows = report_data.read_measurements(
            snapshots["summary"], original_status
        )
        excluded = exclusions(original_status)
        status.update(
            analysis_status=original_status["status"], exclusions=excluded
        )
        save_exclusions(output, excluded)
        for row in excluded:
            log.event(
                "WARNING", "Excluded well", f"{row['Well']}: {row['Reason']}"
            )
        plate = report_data.read_design(
            snapshots["template"],
            args.sheet,
            statistics=args.stats_unit is not None,
            measured_wells=[row["Well"] for row in rows],
        )
        rows, diagnostics = report_data.annotate_wells(rows, plate)
        log.phase(
            2, 4, "Tables and statistics", detail=f"{len(rows)} usable wells"
        )
        warnings = plate["warnings"] + report_data.comparability_notes(rows)
        if original_status.get("error"):
            warnings.append(
                "Cell analysis ended with an error: "
                + original_status["error"]
            )
        no_result = [
            row["Well"] for row in diagnostics if row["Status"] == "NO_RESULT"
        ]
        if no_result:
            warnings.append(
                "Annotated wells without results: " + ", ".join(no_result)
            )
        for warning in warnings:
            log.event("WARNING", "Validation", warning)
        data = {
            **plate,
            "rows": rows,
            "diagnostics": diagnostics,
            "summary": report_data.summarize(rows, plate),
            "warnings": warnings,
            "statistics_enabled": args.stats_unit is not None,
            "comparisons": report_data.comparisons(rows, plate)
            if args.stats_unit
            else [],
            "exclusions": excluded,
            "partial": original_status["status"] == "PARTIAL",
        }
        status.update(
            group_order=plate["group_order"],
            group_order_source=plate["group_order_source"],
            group_order_records=plate["group_order_records"],
        )
        save_json(
            output / "group_order.json",
            {
                key: status[key]
                for key in (
                    "group_order",
                    "group_order_source",
                    "group_order_records",
                )
            },
        )
        log.event(
            "INFO",
            "Group order",
            json.dumps(plate["group_order_records"], ensure_ascii=False),
        )
        for comparison in data["comparisons"]:
            if comparison["Status"] == "Not tested":
                log.event(
                    "WARNING",
                    "Statistics",
                    f"{comparison['Comparison_Block']}: "
                    f"{comparison['Treatment']} versus "
                    f"{comparison['Control']}, {comparison['Metric']}: "
                    f"{comparison['Reason']}",
                )
        log.event(
            "INFO",
            "Statistics",
            "Welch + Holm across both outcomes in each color block"
            if args.stats_unit
            else "Disabled; --stats-unit not supplied",
        )
        save_palette(data, output)
        _save_tables(data, output)
        from .report_output import (
            build_workbook,
            render_plots,
            verify_workbook,
        )

        plots = render_plots(data, output, source.name)
        log.phase(4, 4, "Workbook and verification")
        status["status"] = "PARTIAL" if data["partial"] else "SUCCESS"
        build_workbook(data, plots, snapshots["template"], pending, status)
        verify_workbook(pending, data, plots)
        verify_inputs(manifest)
        final_path = (
            output
            / f"Functional_Report_{safe_label(source.name)}_{run_id}.xlsx"
        )
        pending.replace(final_path)
        plot_records = [
            {
                "Metric": plot["metric"],
                "File": plot["path"].name,
                "Wells": plot["point_wells"],
                "Annotations": plot["annotations"],
                "Palette": "plot_palette.json",
            }
            for plot in plots
        ]
        save_json(output / "plot_manifest.json", plot_records, allow_nan=False)
        status.update(
            status="PARTIAL" if data["partial"] else "SUCCESS",
            workbook=final_path.name,
            measured_wells=len(rows),
            comparison_blocks=len(plate["blocks"]),
            planned_comparisons=len(data["comparisons"]),
            tested_comparisons=sum(
                row["Status"] == "Tested" for row in data["comparisons"]
            ),
            warnings=warnings,
            template_sheet=plate["sheet"],
            inputs=manifest,
            plots=[plot["path"].name for plot in plots],
        )
        log.event(
            status["status"],
            "Report",
            f"{len(rows)} wells; {len(excluded)} excluded; "
            f"workbook: {final_path}",
        )
    except (Exception, KeyboardInterrupt) as error:
        status.update(
            status="CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else ("NO_INPUT" if isinstance(error, NoInputError) else "FAILED"),
            error=str(error) or type(error).__name__,
        )
        pending.unlink(missing_ok=True)
        log.event(status["status"], "Report", status["error"])
        log.event("ERROR", "Traceback", traceback.format_exc(), console=False)
        if isinstance(error, ValidationError) and error.details:
            details = error.details
            columns = list(
                dict.fromkeys(key for row in details for key in row)
            )
            save_csv(
                output / "validation_errors.csv",
                columns,
                details,
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
        folders = read_config(path)
    except (OSError, ValueError) as error:
        command_error("functional_report", error)
        return 1
    succeeded, failed, seen = 0, 0, set()
    for folder in folders:
        canonical = folder.resolve()
        if canonical in seen:
            print(f"Skipping duplicate JSON folder: {folder}")
            continue
        seen.add(canonical)
        try:
            if not folder.is_dir():
                raise NoInputError(f"Folder not found: {folder}")
            result = process_folder(folder, path, args)
            succeeded += result["status"] == "SUCCESS"
            failed += result["status"] != "SUCCESS"
        except KeyboardInterrupt as error:
            command_error("functional_report", error, folder)
            return 130
        except Exception as error:
            failed += 1
            command_error("functional_report", error, folder)
    print(
        f"Functional reports: {succeeded} folder(s) succeeded; "
        f"{failed} partial/failed.",
        flush=True,
    )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
