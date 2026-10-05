"""Join repeated well measurements and compare their baseline changes."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev

from uma_tools.report_schema import ValidationError
from uma_tools.report_statistics import _apply_holm, _welch

from .cell_analysis import SUMMARY_COLUMNS as CELL_COLUMNS
from .report_data import COMPARISON_COLUMNS, METRICS, comparability_notes

RAW_COLUMNS = (
    "Day",
    "Annotation_Status",
    "Comparison_Block",
    "Group",
    "Is_Control",
    "Color_Code",
    "Original_Well",
    *CELL_COLUMNS,
)
CHANGE_COLUMNS = (
    "Day",
    "Baseline_Day",
    "Comparison",
    "Well",
    "Comparison_Block",
    "Group",
    "Is_Control",
    "Metric",
    "Unit",
    "Baseline_Value",
    "Day_Value",
    "Delta",
    "Status",
    "Reason",
)
COVERAGE_COLUMNS = (
    "Day",
    "Well",
    "Excel_Cell",
    "Group",
    "Comparison_Block",
    "Is_Control",
    "Color_Code",
    "Status",
    "Present_In_Baseline",
)
SUMMARY_COLUMNS = (
    "View",
    "Day",
    "Baseline_Day",
    "Comparison_Block",
    "Group",
    "Is_Control",
    "Metric",
    "Unit",
    "Expected_Wells",
    "N",
    "Wells",
    "Mean",
    "SD",
    "Median",
    "Min",
    "Max",
)
STATISTICS_COLUMNS = (
    "Day",
    "Baseline_Day",
    "Comparison",
    "Control_Wells",
    "Treatment_Wells",
    *COMPARISON_COLUMNS,
)
SELECTION_COLUMNS = (
    "Day",
    "Status",
    "Reason",
    "Source_Folder",
    "Selected_Analysis",
    "Summary_File",
    "Completed_Wells",
    "Parameters",
)


def _unique_json(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def _day(value, field):
    if type(value) is not int or value < 0:
        raise ValidationError(f"{field} must be a non-negative integer day")
    return value


def _path(value, parent, field):
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{field} must be a non-empty path string")
    path = Path(value).expanduser()
    return (path if path.is_absolute() else parent / path).resolve()


def read_config(path: Path) -> dict:
    """Resolve paths relative to JSON; never infer days from folder names."""
    if path.name.startswith("._"):
        raise ValidationError("AppleDouble ._ JSON files are not inputs")
    with path.open(encoding="utf-8-sig") as stream:
        config = json.load(stream, object_pairs_hook=_unique_json)
    required = {
        "experiment_name",
        "plate_template",
        "baseline_day",
        "difference_days",
        "timepoints",
    }
    if not isinstance(config, dict):
        raise ValidationError("The survival JSON must contain an object")
    missing = required - config.keys()
    unknown = config.keys() - required - {"output_dir", "sheet"}
    if missing or unknown:
        raise ValidationError(
            f"Invalid survival JSON fields; missing={sorted(missing)}, "
            f"unknown={sorted(unknown)}"
        )
    name = config["experiment_name"]
    if not isinstance(name, str) or not name.strip():
        raise ValidationError("experiment_name must be non-empty text")
    baseline = _day(config["baseline_day"], "baseline_day")
    points = config["timepoints"]
    if not isinstance(points, list) or len(points) < 2:
        raise ValidationError("timepoints must list at least two days")
    timepoints, days, folders = [], set(), set()
    for point in points:
        if not isinstance(point, dict) or set(point) != {"day", "folder"}:
            raise ValidationError(
                "Each timepoint needs exactly day and folder"
            )
        day = _day(point["day"], "timepoints.day")
        folder = _path(point["folder"], path.parent, "timepoints.folder")
        if day in days:
            raise ValidationError(f"Duplicate day: {day}")
        if folder in folders:
            raise ValidationError("Different days cannot use the same folder")
        days.add(day)
        folders.add(folder)
        timepoints.append({"day": day, "folder": folder})
    if baseline not in days:
        raise ValidationError("baseline_day is not present in timepoints")
    differences = config["difference_days"]
    if not isinstance(differences, list) or not differences:
        raise ValidationError("difference_days must be a non-empty list")
    differences = [_day(day, "difference_days") for day in differences]
    if len(set(differences)) != len(differences):
        raise ValidationError("Duplicate difference_days")
    if baseline in differences or not set(differences) <= days:
        raise ValidationError(
            "difference_days must refer to recorded days other than baseline"
        )
    template = _path(config["plate_template"], path.parent, "plate_template")
    if template.suffix.lower() != ".xlsx" or template.name.startswith(
        (".", "~$")
    ):
        raise ValidationError("plate_template must be a regular .xlsx input")
    sheet = config.get("sheet")
    if sheet is not None and (not isinstance(sheet, str) or not sheet.strip()):
        raise ValidationError("sheet must be a non-empty worksheet name")
    return {
        **config,
        "experiment_name": name.strip(),
        "plate_template": template,
        "output_dir": _path(
            config.get("output_dir", "."), path.parent, "output_dir"
        ),
        "timepoints": sorted(timepoints, key=lambda point: point["day"]),
        "difference_days": sorted(differences),
        "sheet": sheet,
    }


def join_days(measurements: dict[int, list], plate: dict, baseline: int):
    """Retain all input rows, including wells outside the annotated design."""
    design = {row["Well"]: row for row in plate["design"]}
    baseline_wells = {row["Well"] for row in measurements[baseline]}
    records, coverage, warnings = [], [], []
    for day, rows in sorted(measurements.items()):
        measured = {row["Well"] for row in rows}
        unmapped = sorted(measured - design.keys())
        missing = sorted(design.keys() - measured)
        additional = sorted(measured - baseline_wells)
        if additional:
            warnings.append(
                f"Day {day}: additional wells relative to baseline Day "
                f"{baseline}: {', '.join(additional)}"
            )
        if unmapped:
            warnings.append(
                f"Day {day}: unannotated wells retained in raw data, excluded "
                f"from grouped plots and statistics: {', '.join(unmapped)}"
            )
        if missing:
            warnings.append(
                f"Day {day}: annotated wells without measurements: "
                + ", ".join(missing)
            )
        for row in rows:
            annotation = design.get(row["Well"], {})
            records.append(
                {
                    **row,
                    "Day": day,
                    "Annotation_Status": "MAPPED"
                    if annotation
                    else "UNMAPPED",
                    **{
                        key: annotation.get(key)
                        for key in (
                            "Comparison_Block",
                            "Group",
                            "Is_Control",
                            "Color_Code",
                        )
                    },
                }
            )
        for well, coordinate in plate["cells"].items():
            annotation = design.get(well, {})
            state = (
                ("MEASURED" if annotation else "UNMAPPED")
                if well in measured
                else ("NO_RESULT" if annotation else "UNUSED")
            )
            coverage.append(
                {
                    "Day": day,
                    "Well": well,
                    "Excel_Cell": f"{plate['sheet']}!{coordinate}",
                    "Status": state,
                    "Present_In_Baseline": well in baseline_wells,
                    **{
                        key: annotation.get(key)
                        for key in (
                            "Comparison_Block",
                            "Group",
                            "Is_Control",
                            "Color_Code",
                        )
                    },
                }
            )
    warnings.extend(comparability_notes(records))
    methods = {row["Threshold_Method"] for row in records}
    if len(methods) > 1:
        warnings.append("Threshold methods differ across days/wells")
    manual = {
        (row["Threshold_Lower"], row["Threshold_Upper"])
        for row in records
        if row["Threshold_Method"] == "manual"
    }
    if len(manual) > 1:
        warnings.append("Manual threshold bounds differ across days/wells")
    return records, coverage, warnings


def calculate_changes(rows, plate, baseline, difference_days):
    """Subtract matched wells; never replace missing values with zero."""
    lookup = {(row["Day"], row["Well"]): row for row in rows}
    if len(lookup) != len(rows):
        raise ValidationError("Duplicate day/well measurements")
    changes = []
    for day in difference_days:
        for annotation in plate["design"]:
            well = annotation["Well"]
            before, after = (
                lookup.get((baseline, well)),
                lookup.get((day, well)),
            )
            if before is not None and after is not None:
                state, reason = "PAIRED", ""
            elif before is None and after is None:
                state, reason = "MISSING_BOTH", "No measurements in either day"
            elif before is None:
                state, reason = "MISSING_BASELINE", "No baseline measurement"
            else:
                state, reason = "MISSING_DAY", "No comparison-day measurement"
            for metric, _, unit in METRICS:
                a = before[metric] if before is not None else None
                b = after[metric] if after is not None else None
                changes.append(
                    {
                        "Day": day,
                        "Baseline_Day": baseline,
                        "Comparison": f"{day}-{baseline}",
                        "Well": well,
                        **{
                            key: annotation[key]
                            for key in (
                                "Comparison_Block",
                                "Group",
                                "Is_Control",
                            )
                        },
                        "Metric": metric,
                        "Unit": unit,
                        "Baseline_Value": a,
                        "Day_Value": b,
                        "Delta": b - a if state == "PAIRED" else None,
                        "Status": state,
                        "Reason": reason,
                    }
                )
    return changes


def observations(data, view, metric):
    """Return only eligible, measured well values for a plot."""
    if view == "changes":
        return [
            {**row, "Value": row["Delta"]}
            for row in data["changes"]
            if row["Metric"] == metric and row["Status"] == "PAIRED"
        ]
    if view not in {"raw", "baseline"}:
        raise ValueError(f"Unknown plot view: {view}")
    return [
        {**row, "Value": row[metric]}
        for row in data["rows"]
        if row["Annotation_Status"] == "MAPPED"
        and (view == "raw" or row["Day"] == data["baseline_day"])
    ]


def summarize(data):
    result = []
    for view, days in (
        ("raw", data["days"]),
        ("changes", data["difference_days"]),
    ):
        for metric, _, unit in METRICS:
            observed = observations(data, view, metric)
            for block in data["blocks"]:
                for group in block["groups"]:
                    expected = [
                        row
                        for row in data["design"]
                        if row["Comparison_Block"] == block["id"]
                        and row["Group"] == group
                    ]
                    for day in days:
                        selected = [
                            row
                            for row in observed
                            if row["Comparison_Block"] == block["id"]
                            and row["Group"] == group
                            and row["Day"] == day
                        ]
                        values = [row["Value"] for row in selected]
                        result.append(
                            {
                                "View": view,
                                "Day": day,
                                "Baseline_Day": data["baseline_day"]
                                if view == "changes"
                                else None,
                                "Comparison_Block": block["id"],
                                "Group": group,
                                "Is_Control": expected[0]["Is_Control"],
                                "Metric": metric,
                                "Unit": unit,
                                "Expected_Wells": len(expected),
                                "N": len(values),
                                "Wells": ", ".join(
                                    row["Well"] for row in selected
                                ),
                                "Mean": mean(values) if values else None,
                                "SD": stdev(values)
                                if len(values) >= 2
                                else None,
                                "Median": median(values) if values else None,
                                "Min": min(values) if values else None,
                                "Max": max(values) if values else None,
                            }
                        )
    return result


def compare_changes(data):
    """Welch on per-well changes; one planned Holm family per color."""
    grouped = defaultdict(list)
    for row in data["changes"]:
        if row["Status"] == "PAIRED":
            grouped[
                (
                    row["Comparison_Block"],
                    row["Group"],
                    row["Day"],
                    row["Metric"],
                )
            ].append(row)
    comparisons = []
    for block in data["blocks"]:
        control = block["control"]
        family_size = (
            len(METRICS)
            * len(data["difference_days"])
            * (len(block["groups"]) - 1)
        )
        for day in data["difference_days"]:
            for group in block["groups"]:
                if group == control:
                    continue
                for metric, _, unit in METRICS:
                    treatments = grouped[(block["id"], group, day, metric)]
                    controls = grouped[(block["id"], control, day, metric)]
                    a = [row["Delta"] for row in treatments]
                    b = [row["Delta"] for row in controls]
                    test, reason = _welch(a, b)
                    result = dict.fromkeys(STATISTICS_COLUMNS)
                    result.update(
                        {
                            "Day": day,
                            "Baseline_Day": data["baseline_day"],
                            "Comparison": f"{day}-{data['baseline_day']}",
                            "Comparison_Block": block["id"],
                            "Color_Code": block["color_key"],
                            "Control": control,
                            "Treatment": group,
                            "Metric": metric,
                            "Unit": unit,
                            "Stats_Unit": "well",
                            "Control_N": len(b),
                            "Treatment_N": len(a),
                            "Control_Wells": ", ".join(
                                r["Well"] for r in controls
                            ),
                            "Treatment_Wells": ", ".join(
                                r["Well"] for r in treatments
                            ),
                            "Control_Mean": mean(b) if b else None,
                            "Treatment_Mean": mean(a) if a else None,
                            "Difference": mean(a) - mean(b)
                            if a and b
                            else None,
                            "Family_Size": family_size,
                            "Status": "Tested" if test else "Not tested",
                            "Significance": "",
                            "Reason": reason,
                            **test,
                        }
                    )
                    comparisons.append(result)
    _apply_holm(comparisons)
    return comparisons
