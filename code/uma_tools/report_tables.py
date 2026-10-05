"""
Read literal summary tables, source-image identities, and plate maps.
"""

from __future__ import annotations

import csv
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

from .report_schema import (
    ALIGNMENT_FILENAME_SUFFIX,
    FN_METRIC,
    NUMBER_PATTERN,
    POINT_PATTERN,
    POINT_WELL_PATTERN,
    SEQUENCE_PATTERN,
    WELL_PATTERN,
    ParsedRecord,
    ValidationError,
)


def read_csv_table(
    path: Path, label: str
) -> tuple[list[str], list[ParsedRecord]]:
    """
    Read comma-delimited UTF-8 CSV without guessing delimiters or
    skipping rows.
    """
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, strict=True)
        columns = next(reader, None)
        if not columns or any(not value.strip() for value in columns):
            raise ValidationError(f"{label}: missing or blank CSV headers.")
        if len(columns) != len(set(columns)):
            raise ValidationError(f"{label}: duplicate CSV headers.")
        if "File_Name" not in columns:
            raise ValidationError(
                f"{label}: required column File_Name was not found. "
                "Expected comma-delimited CSV."
            )
        records = []
        for source_row, fields in enumerate(reader, 2):
            if len(fields) != len(columns):
                raise ValidationError(
                    f"{label}: invalid field count in CSV "
                    f"record {source_row}.",
                    [
                        {
                            "Table": label,
                            "Source_Row": source_row,
                            "Issue": (
                                f"Expected {len(columns)} fields; "
                                f"found {len(fields)}"
                            ),
                        }
                    ],
                )
            records.append(
                {"source_row": source_row, "raw": dict(zip(columns, fields))}
            )
    if not records:
        raise ValidationError(f"{label}: the table contains no image records.")
    return columns, records


def parse_filenames(
    records: list[ParsedRecord],
    label: str,
    original_stems: dict[str, str] | None = None,
) -> None:
    """
    Resolve exact identities and parse optional acquisition metadata.

    Original filenames are literal identifiers, not paths. Only
    alignment's documented terminal suffix is removed to look up the
    original complete stem. A Seq token never changes identity. Dots,
    Unicode, and terminal Well tokens are accepted; Seq and Point
    metadata are optional.
    """
    errors = []
    for record in records:
        name = record["raw"]["File_Name"]
        issue = None
        image_id = name
        if not name or not name.strip() or "/" in name or "\\" in name:
            issue = (
                "File_Name must be a nonempty original filename, not a path."
            )
        elif original_stems is not None:
            if not name.endswith(ALIGNMENT_FILENAME_SUFFIX):
                issue = (
                    "Alignment File_Name must end with "
                    f"{ALIGNMENT_FILENAME_SUFFIX!r}."
                )
            else:
                stem = name[: -len(ALIGNMENT_FILENAME_SUFFIX)]
                image_id = original_stems.get(stem)
                if image_id is None:
                    issue = (
                        "Alignment source stem has no exact original "
                        "File_Name match in thickness/area."
                    )
        if issue is None:
            original_stem = Path(image_id).stem
            wells = list(WELL_PATTERN.finditer(original_stem))
            points = list(POINT_WELL_PATTERN.finditer(original_stem))
            numbers = list(POINT_PATTERN.finditer(original_stem))
            sequences = list(SEQUENCE_PATTERN.finditer(original_stem))
            if len(wells) != 1:
                issue = (
                    "Expected exactly one valid 96-well WellA1/WellA01 token."
                )
            elif len(points) > 1 or len(sequences) > 1:
                issue = "Repeated Point or Seq metadata tokens are ambiguous."
            else:
                well = f"{wells[0][1].upper()}{int(wells[0][2]):02d}"
                point_well = (
                    f"{points[0][1].upper()}{int(points[0][2]):02d}"
                    if points
                    else None
                )
                if point_well is not None and well != point_well:
                    issue = f"Well ({well}) and Point ({point_well}) disagree."
                else:
                    record.update(
                        image_id=image_id,
                        well=well,
                        image_number=numbers[0][3].zfill(4)
                        if numbers
                        else None,
                        sequence_number=sequences[0][1] if sequences else None,
                    )
        if issue:
            errors.append(
                {
                    "Table": label,
                    "Source_Row": record["source_row"],
                    "File_Name": name,
                    "Issue": issue,
                }
            )
    if errors:
        raise ValidationError(
            f"{label}: {len(errors)} filename(s) could not be parsed.", errors
        )
    counts = Counter(record["image_id"] for record in records)
    duplicates = [
        {
            "Table": label,
            "Source_Row": row["source_row"],
            "Image_ID": row["image_id"],
            "File_Name": row["raw"]["File_Name"],
            "Issue": "Duplicate Image_ID",
        }
        for row in records
        if counts[row["image_id"]] > 1
    ]
    if duplicates:
        raise ValidationError(
            f"{label}: duplicate Image_ID values. No rows were removed.",
            duplicates,
        )


def original_stem_lookup(thickness, fibronectin):
    """
    Reject extension collisions before resolving alignment source stems.
    """
    by_stem = {}
    for record in thickness + fibronectin:
        image_id = record["image_id"]
        by_stem.setdefault(Path(image_id).stem, set()).add(image_id)
    ambiguous = [
        {
            "Source_Stem": stem,
            "Image_ID": image_id,
            "Issue": (
                "Ambiguous original stem: alignment filenames omit the source "
                "extension"
            ),
        }
        for stem, image_ids in sorted(by_stem.items())
        if len(image_ids) > 1
        for image_id in sorted(image_ids)
    ]
    if ambiguous:
        raise ValidationError(
            (
                "Original filenames share a stem; alignment cannot identify "
                "their extensions."
            ),
            ambiguous,
        )
    return {stem: next(iter(image_ids)) for stem, image_ids in by_stem.items()}


def optional_number_sort(value):
    """
    Order optional metadata numerically, with absent values after
    numbers.
    """
    return (value is None, 0 if value is None else int(value))


def numeric_values(records, fields, label, upper=None):
    errors = []
    for record in records:
        record["numbers"] = {}
        for field in fields:
            original = record["raw"][field]
            value = (
                float(original)
                if NUMBER_PATTERN.fullmatch(original.strip())
                else math.nan
            )
            if (
                not math.isfinite(value)
                or value < 0
                or (upper is not None and value > upper)
            ):
                errors.append(
                    {
                        "Table": label,
                        "Source_Row": record["source_row"],
                        "Image_ID": record["image_id"],
                        "Metric": field,
                        "Original_Value": original,
                        "Issue": (
                            "Missing, non-numeric, non-finite, negative, or "
                            "out-of-range value"
                        ),
                    }
                )
            else:
                record["numbers"][field] = value
    if errors:
        raise ValidationError(
            f"{label}: invalid required measurements. "
            "No values were imputed or excluded.",
            errors,
        )


def metadata_column(values):
    """
    Keep text and identifiers intact; type entirely numeric metadata
    columns.
    """
    populated = [value for value in values if value != ""]
    is_number = bool(populated) and all(
        NUMBER_PATTERN.fullmatch(value.strip()) for value in populated
    )
    if is_number and not any(
        re.match(r"^[+-]?0[0-9]", value.strip()) for value in populated
    ):
        numbers = [None if value == "" else float(value) for value in values]
        if all(
            value is None or (math.isfinite(value) and abs(value) < 10**15)
            for value in numbers
        ):
            return [
                int(value)
                if value is not None and value.is_integer()
                else value
                for value in numbers
            ]
    return [None if value == "" else value for value in values]


def validate_fn_threshold(value: object) -> float:
    """
    The explicit percentage cutoff is inclusive for retained
    observations.
    """
    try:
        threshold = float(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(
            "FN_AREA_THRESHOLD_PERCENT must be a number from 0 to 100."
        ) from error
    if (
        isinstance(value, bool)
        or not math.isfinite(threshold)
        or not 0 <= threshold <= 100
    ):
        raise ValidationError(
            (
                "FN_AREA_THRESHOLD_PERCENT must be finite and between 0 and "
                "100 inclusive."
            )
        )
    return threshold


def validate_fibronectin(records, columns):
    """
    Validate reported coverage and reconcile optional source
    identifiers/counts.
    """
    if FN_METRIC not in columns:
        raise ValidationError(
            f"Fibronectin is missing required column: {FN_METRIC}"
        )
    numeric_values(records, [FN_METRIC], "Fibronectin", upper=100)
    errors = []
    reconcile_pixels = {"FN_Positive_Pixels", "Total_Pixels"}.issubset(columns)
    for row in records:
        raw, image_id = row["raw"], row["image_id"]
        if "Image_ID" in columns and raw["Image_ID"] not in (
            image_id,
            Path(image_id).stem,
        ):
            errors.append(
                {
                    "Table": "Fibronectin",
                    "Source_Row": row["source_row"],
                    "Image_ID": image_id,
                    "Source_Image_ID": raw["Image_ID"],
                    "Issue": (
                        "Image_ID must equal the complete File_Name or its "
                        "exact legacy source stem"
                    ),
                }
            )
        if reconcile_pixels:
            try:
                positive, total = (
                    float(raw["FN_Positive_Pixels"]),
                    float(raw["Total_Pixels"]),
                )
                valid = (
                    math.isfinite(positive)
                    and math.isfinite(total)
                    and positive.is_integer()
                    and total.is_integer()
                    and 0 <= positive <= total
                    and total > 0
                )
                valid = valid and math.isclose(
                    row["numbers"][FN_METRIC],
                    positive / total * 100,
                    rel_tol=1e-9,
                    abs_tol=1e-6,
                )
            except (ValueError, TypeError, OverflowError):
                valid = False
            if not valid:
                errors.append(
                    {
                        "Table": "Fibronectin",
                        "Source_Row": row["source_row"],
                        "Image_ID": image_id,
                        "Issue": (
                            "FN_Area_Percent or pixel counts are inconsistent "
                            "with 100 * positive / total"
                        ),
                    }
                )
    if errors:
        raise ValidationError(
            "Fibronectin source identifiers or pixel counts are inconsistent.",
            errors,
        )
    return reconcile_pixels


def _mask_bound(record, column):
    """Read an optional effective bound without inventing defaults."""
    value = record["raw"].get(column, "").strip()
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        number = math.nan
    if not math.isfinite(number) or number < 0:
        raise ValidationError(
            f"Invalid FN mask {column} in source row {record['source_row']}: "
            f"{value!r}",
            [
                {
                    "Table": "Fibronectin",
                    "Source_Row": record["source_row"],
                    "File_Name": record["raw"]["File_Name"],
                    "Issue": f"Invalid {column}: {value!r}",
                }
            ],
        )
    return number


def _mask_description(projection, lower, upper, units):
    """Describe source mask bounds separately from a coverage cutoff."""

    def bound(value):
        if value is None:
            return "not recorded"
        if value == float.fromhex("0x1.fffffep+127"):
            return "float32 max"
        return f"{value:.9g}"

    method = (
        f"{projection} projection" if projection else "projection not recorded"
    )
    if lower is None and upper is None:
        return f"{method}; intensity thresholds not recorded"
    intensity = (
        "raw intensity"
        if units.lower() == "raw projection intensity"
        else units or "intensity (units not recorded)"
    )
    if lower is None or upper is None:
        return (
            f"{method}; {intensity} lower={bound(lower)}, upper={bound(upper)}"
        )
    return (
        f"{method}; {intensity} [{bound(lower)}, {bound(upper)}] (inclusive)"
    )


def summarize_fn_mask_settings(records, log):
    """Label FN masks from area CSV metadata, not report defaults."""
    profiles = {}
    for record in records:
        raw = record["raw"]
        lower = _mask_bound(record, "Threshold_Lower")
        upper_column = (
            "Effective_Threshold_Upper"
            if raw.get("Effective_Threshold_Upper", "").strip()
            else "Threshold_Upper"
        )
        upper = _mask_bound(record, upper_column)
        if lower is not None and upper is not None and upper < lower:
            raise ValidationError(
                "FN mask upper threshold is below its lower threshold in "
                f"source row {record['source_row']} ({raw['File_Name']})."
            )
        projection = raw.get("Projection_Method", "").strip().upper()
        units = raw.get("Threshold_Units", "").strip()
        key = (projection, lower, upper, units)
        if key not in profiles:
            profiles[key] = {
                "projection": projection or None,
                "lower": lower,
                "upper": upper,
                "units": units or None,
                "images": 0,
                "description": _mask_description(*key),
            }
        profiles[key]["images"] += 1
    values = list(profiles.values())
    mixed = len(values) > 1
    complete = all(
        item["lower"] is not None
        and item["upper"] is not None
        and item["projection"]
        and item["units"]
        for item in values
    )
    caption = (
        f"FN mask: mixed intensity/projection settings ({len(values)}); "
        "see Overview and FN_Source columns."
        if mixed
        else f"FN mask: {values[0]['description']}."
    )
    log.event(
        "PASS" if complete and not mixed else "WARNING",
        "FN mask thresholds",
        caption,
    )
    return {
        "caption": caption,
        "profiles": values,
        "complete": complete,
        "mixed": mixed,
    }


def read_template(
    path: Path, sheet_name: str | None
) -> tuple[str, list[list[Any]], dict[str, str], dict[str, str]]:
    """
    Read literal group labels at their real Excel coordinates, without
    shifting.
    """
    import openpyxl

    workbook = openpyxl.load_workbook(path, read_only=False, data_only=False)
    try:
        selected = workbook.sheetnames[0] if sheet_name is None else sheet_name
        if selected not in workbook.sheetnames:
            raise ValidationError(
                f"Template worksheet {selected!r} was not found. "
                f"Available: {workbook.sheetnames}"
            )
        sheet = workbook[selected]
        for merged in sheet.merged_cells.ranges:
            if merged.min_row <= 9 and merged.min_col <= 13:
                raise ValidationError(
                    f"Template contains merged cells within A1:M9: {merged}. "
                    "Use one group per well cell."
                )
        grid = [
            [sheet.cell(row, column).value for column in range(1, 14)]
            for row in range(1, 10)
        ]
        for index, value in enumerate(grid[0][1:], 1):
            if isinstance(value, bool) or str(value).strip() not in (
                str(index),
                f"{index}.0",
            ):
                raise ValidationError(
                    (
                        "Invalid template headers. B1:M1 must contain columns "
                        "1-12 in order."
                    )
                )
        for index, letter in enumerate("ABCDEFGH", 1):
            if str(grid[index][0]).strip().upper() != letter:
                raise ValidationError(
                    (
                        "Invalid template headers. A2:A9 must contain rows "
                        "A-H in order."
                    )
                )
        well_map, cells = {}, {}
        for row_index, letter in enumerate("ABCDEFGH", 2):
            for column in range(1, 13):
                well = f"{letter}{column:02d}"
                cell = sheet.cell(row_index, column + 1)
                cells[well] = cell.coordinate
                if cell.data_type in ("f", "e"):
                    raise ValidationError(
                        f"Template {selected}!{cell.coordinate} must contain "
                        "a literal group name, not a formula or Excel error."
                    )
                value = cell.value
                if value is not None and str(value).strip():
                    well_map[well] = str(value)
        if not well_map:
            raise ValidationError("The template contains no annotated wells.")
        return selected, grid, well_map, cells
    finally:
        workbook.close()


def prepare_display_tables(data, plots, parameters, manifests):
    """Expose plot observations and provenance without recalculation."""
    import json

    tables = {}

    def add(name, columns, rows):
        tables[name] = {"columns": columns, "rows": rows}

    overview = [
        ("Report", "UMA fibronectin analysis"),
        ("Source folder", parameters.get("source_folder", "")),
        ("Collection", parameters.get("combined_results_folder", "")),
        ("Plate", data["plate_id"]),
        ("Run ID", parameters.get("run_id", "")),
        ("Images before FN filtering", len(data["rows"])),
        ("Images retained", len(data["retained_rows"])),
        ("Images below FN cutoff", len(data["excluded_rows"])),
        (
            "Registered processing exclusions",
            len(data.get("processing_exclusions", [])),
        ),
        ("FN cutoff (%)", data["fn_threshold"]),
        ("FN mask intensity thresholds", data["fn_mask_settings"]["caption"]),
        ("Statistics unit", parameters.get("stats_unit") or "Disabled"),
        ("Figures", "Seven full-data and seven filtered views on Plots"),
        ("Point color", "Condition; lavender reserved for a bold control"),
        ("Point shape", "Technical well within each condition; see Plot_Data"),
        ("Plot format", parameters.get("plot_format", "pdf")),
        (
            "Population",
            "Points and boxes represent images; "
            "wells are technical replicates.",
        ),
        (
            "Snapshot",
            "Rerun the report after changing inputs; "
            "embedded plots do not recalculate.",
        ),
    ]
    overview.extend(
        (
            f"FN mask setting {index}",
            f"{profile['description']}; {profile['images']} source image(s)",
        )
        for index, profile in enumerate(
            data["fn_mask_settings"]["profiles"], 1
        )
    )
    add(
        "Overview",
        ["Item", "Value"],
        [dict(zip(("Item", "Value"), r)) for r in overview],
    )
    observations, labels, info = [], [], []
    for plot in plots:
        row_lookup = {row["Image_ID"]: row for row in data["rows"]}
        memberships = {
            g: panel["id"] for panel in plot["panels"] for g in panel["groups"]
        }
        for image_id in plot["plotted_image_ids"]:
            row = row_lookup[image_id]
            replicate = data["group_wells"][row["Group"]].index(row["Well"])
            observations.append(
                {
                    "Plot": plot["plot_id"],
                    "View": plot["view"],
                    "Metric": plot["metric"],
                    "Unit": plot["unit"],
                    "Panel": memberships[row["Group"]],
                    "Group": row["Group"],
                    "Well": row["Well"],
                    "Image_ID": image_id,
                    "Value": row[plot["metric"]],
                    "Technical_Replicate": replicate + 1,
                    "Point_Color": plot["condition_styles"][row["Group"]][
                        "color"
                    ],
                    "Point_Marker": plot["well_markers"][row["Well"]],
                    "Red_Outline": image_id in plot["red_outline_image_ids"],
                    "Panel_X": plot["rendered_x_positions"][image_id],
                }
            )
        for panel in plot["panels"]:
            for group in panel["groups"]:
                labels.append(
                    {
                        "Plot": plot["plot_id"],
                        "Panel": panel["id"],
                        "Color_Code": panel["color_code"],
                        "Figure_Context": panel["context"],
                        "Panel_Title": panel["title"],
                        "Group": group,
                        "Display_Label": panel["labels"][group],
                        "Point_Color": plot["condition_styles"][group][
                            "color"
                        ],
                        "Color_Name": plot["condition_styles"][group][
                            "color_name"
                        ],
                        "Is_Control": plot["condition_styles"][group][
                            "is_control"
                        ],
                        "Palette_Mode": plot["condition_styles"][group][
                            "palette_mode"
                        ],
                        "N_Images": plot["group_counts"][group],
                        "N_Wells": plot["group_well_counts"][group],
                    }
                )
        info.append(
            {
                "Plot": plot["plot_id"],
                "Metric": plot["metric"],
                "View": plot["view"],
                "Title": plot["title"],
                "File": Path(plot["path"]).name,
                "PDF_File": Path(plot["pdf_file"]).name
                if plot["pdf_file"]
                else "",
                "PNG_File": Path(plot["png_file"]).name
                if plot["png_file"]
                else "",
                "Font": plot["font"],
                "PNG_DPI": plot["png_dpi"],
                "Panels": len(plot["panels"]),
                "Images": plot["point_count"],
                "Y_Min": plot["y_min"],
                "Y_Max": plot["y_max"],
                "Caption": plot["caption"],
                "Statistics": plot["statistics_note"],
            }
        )
    add(
        "Plot_Data",
        [
            "Plot",
            "View",
            "Metric",
            "Unit",
            "Panel",
            "Group",
            "Well",
            "Image_ID",
            "Value",
            "Technical_Replicate",
            "Point_Color",
            "Point_Marker",
            "Red_Outline",
            "Panel_X",
        ],
        observations,
    )
    add(
        "Plot_Labels",
        [
            "Plot",
            "Panel",
            "Color_Code",
            "Figure_Context",
            "Panel_Title",
            "Group",
            "Display_Label",
            "Point_Color",
            "Color_Name",
            "Is_Control",
            "Palette_Mode",
            "N_Images",
            "N_Wells",
        ],
        labels,
    )
    add(
        "Plot_Info",
        [
            "Plot",
            "Metric",
            "View",
            "Title",
            "File",
            "PDF_File",
            "PNG_File",
            "Font",
            "PNG_DPI",
            "Panels",
            "Images",
            "Y_Min",
            "Y_Max",
            "Caption",
            "Statistics",
        ],
        info,
    )
    run_rows = [
        {
            "Parameter": key,
            "Value": json.dumps(value, ensure_ascii=False)
            if isinstance(value, (dict, list, tuple))
            else value,
        }
        for key, value in parameters.items()
        if key not in ("status", "stage")
    ]
    add("Run_Info", ["Parameter", "Value"], run_rows)
    add(
        "Source_Files",
        [
            "Input",
            "Path",
            "Archived_Path",
            "Archive_Relative_Path",
            "Bytes",
            "SHA256",
        ],
        manifests,
    )
    data["display_tables"] = tables
