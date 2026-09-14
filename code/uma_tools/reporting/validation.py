"""
Reconcile all assays, apply plate annotations, and derive FN-filtered
views.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from .constants import (
    ALIGNMENT_PATTERN,
    FN_INCLUDED_FLAG,
    FN_LOW_FLAG,
    FN_METRIC,
    FN_REASON_COLUMN,
    FN_THRESHOLD_COLUMN,
    THICKNESS_METRICS,
    THICKNESS_UNITS,
)
from .io import save_csv
from .models import EventLogger, ReportData, ReportInputs, ValidationError
from .plate import read_template
from .source_tables import (
    metadata_column,
    numeric_values,
    optional_number_sort,
    original_stem_lookup,
    parse_filenames,
    read_csv_table,
    validate_fibronectin,
    validate_fn_threshold,
)


def filter_counts(rows, group_order, group_wells):
    """
    Count image-level exclusions while keeping the original group/well
    order.
    """
    retained = [row for row in rows if row[FN_INCLUDED_FLAG]]
    excluded = [row for row in rows if row[FN_LOW_FLAG]]
    if len(retained) + len(excluded) != len(rows):
        raise RuntimeError(
            "Full, retained, and excluded image counts do not reconcile."
        )
    groups, wells = [], []
    for group in group_order:
        group_rows = [row for row in rows if row["Group"] == group]
        kept = [row for row in group_rows if row[FN_INCLUDED_FLAG]]
        groups.append(
            {
                "Group": group,
                "Total_Images": len(group_rows),
                "Retained_Images": len(kept),
                "Excluded_Images": len(group_rows) - len(kept),
                "Original_Wells": len(group_wells[group]),
                "Wells_With_Retained_Images": len(
                    {row["Well"] for row in kept}
                ),
            }
        )
        for replicate, well in enumerate(group_wells[group], 1):
            original_count = sum(row["Well"] == well for row in group_rows)
            retained_count = sum(row["Well"] == well for row in kept)
            wells.append(
                {
                    "Group": group,
                    "Well": well,
                    "Technical_Replicate": replicate,
                    "Total_Images": original_count,
                    "Retained_Images": retained_count,
                    "Excluded_Images": original_count - retained_count,
                }
            )
    return retained, excluded, groups, wells


def validate_and_merge(
    paths: ReportInputs,
    sheet_name: str | None,
    plate_label: str,
    output: Path,
    log: EventLogger,
    fn_threshold: float,
) -> ReportData:
    fn_threshold = validate_fn_threshold(fn_threshold)
    alignment_columns, alignment = read_csv_table(
        paths["alignment"], "Alignment"
    )
    thickness_columns, thickness = read_csv_table(
        paths["thickness"], "Thickness"
    )
    fn_columns, fibronectin = read_csv_table(
        paths["fibronectin"], "Fibronectin"
    )
    log.event(
        "INFO",
        "Input rows",
        f"Alignment: {len(alignment)}; thickness: {len(thickness)}; "
        f"fibronectin: {len(fibronectin)}",
    )
    matches = [
        (column, ALIGNMENT_PATTERN.fullmatch(column))
        for column in alignment_columns
        if ALIGNMENT_PATTERN.fullmatch(column)
    ]
    if len(matches) != 1:
        raise ValidationError(
            "Expected exactly one alignment-percentage column; "
            f"found {len(matches)}."
        )
    metric, match = matches[0]
    angle_label = match[1]
    missing = [
        field for field in THICKNESS_METRICS if field not in thickness_columns
    ]
    if missing:
        raise ValidationError(
            "Thickness is missing required columns: " + ", ".join(missing)
        )
    parse_filenames(thickness, "Thickness")
    parse_filenames(fibronectin, "Fibronectin")
    parse_filenames(
        alignment, "Alignment", original_stem_lookup(thickness, fibronectin)
    )
    alignment_ids = {row["image_id"] for row in alignment}
    thickness_ids = {row["image_id"] for row in thickness}
    fn_ids = {row["image_id"] for row in fibronectin}
    if not alignment_ids == thickness_ids == fn_ids:
        union = alignment_ids | thickness_ids | fn_ids
        details = [
            {
                "Table": source,
                "Image_ID": value,
                "Issue": "Image_ID is absent from this table",
            }
            for source, values in (
                ("Alignment", alignment_ids),
                ("Thickness", thickness_ids),
                ("Fibronectin", fn_ids),
            )
            for value in sorted(union - values)
        ]
        raise ValidationError(
            (
                "The three CSV tables do not have identical one-to-one "
                "Image_ID sets."
            ),
            details,
        )
    log.event(
        "PASS",
        "Image matching",
        f"{len(alignment_ids)}/{len(alignment_ids)} unique images matched "
        "across all three CSV tables.",
    )
    selected, grid, well_map, cells = read_template(
        paths["template"], sheet_name
    )
    log.event(
        "INFO", "Template", f"{paths['template']} | Worksheet: {selected}"
    )
    log.event("INFO", "Annotated wells", ", ".join(well_map))
    image_counts = Counter(row["well"] for row in alignment)
    diagnostic = [
        {
            "Well": well,
            "Template_Cell": f"{selected}!{cell}",
            "Group": well_map.get(well, ""),
            "Image_Count": image_counts[well],
            "Status": ("UNANNOTATED" if well not in well_map else "MAPPED")
            if image_counts[well]
            else ("NO_IMAGES" if well in well_map else "EMPTY"),
        }
        for well, cell in cells.items()
    ]
    save_csv(
        output / "annotation_diagnostics.csv", list(diagnostic[0]), diagnostic
    )
    uncovered = [row for row in alignment if row["well"] not in well_map]
    matched_count = len(alignment) - len(uncovered)
    coverage = 100 * matched_count / len(alignment)
    log.event(
        "INFO",
        "Annotation coverage",
        f"{matched_count}/{len(alignment)} images ({coverage:.2f}%)",
    )
    if uncovered:
        wells = sorted({row["well"] for row in uncovered})
        details = [
            {
                "Table": "Alignment",
                "Source_Row": row["source_row"],
                "Image_ID": row["image_id"],
                "Well": row["well"],
                "Template_Cell": f"{selected}!{cells[row['well']]}",
                "Issue": "No group assigned in the selected template",
            }
            for row in uncovered
        ]
        raise ValidationError(
            f"Annotation coverage is {coverage:.2f}% "
            f"({matched_count}/{len(alignment)}), not 100%. "
            f"Template: {paths['template'].name}. "
            f"Unannotated wells: {', '.join(wells)}. "
            "Check the selected template workbook and the indicated Excel "
            "cells. No automatic template shifting is performed.",
            details,
        )
    unused_wells = [well for well in well_map if not image_counts[well]]
    if unused_wells:
        log.event(
            "WARNING", "Template wells without images", ", ".join(unused_wells)
        )
    numeric_values(alignment, [metric], "Alignment", upper=100)
    numeric_values(thickness, THICKNESS_METRICS, "Thickness")
    pixel_counts_checked = validate_fibronectin(fibronectin, fn_columns)
    inconsistent = [
        {
            "Table": "Thickness",
            "Source_Row": row["source_row"],
            "Image_ID": row["image_id"],
            **row["numbers"],
            "Issue": "Min <= Median <= Max is not satisfied",
        }
        for row in thickness
        if not row["numbers"]["Min"]
        <= row["numbers"]["Median"]
        <= row["numbers"]["Max"]
    ]
    if inconsistent:
        raise ValidationError(
            "Thickness contains inconsistent Min, Median, and Max values.",
            inconsistent,
        )
    log.event(
        "PASS",
        "Numeric validation",
        (
            "All seven metrics are complete, finite, and within the required "
            "ranges."
        ),
    )

    used_groups = {well_map[well] for well in image_counts}
    group_order = [
        group
        for group in dict.fromkeys(well_map.values())
        if group in used_groups
    ]
    group_wells = {
        group: [
            well
            for well in well_map
            if image_counts[well] and well_map[well] == group
        ]
        for group in group_order
    }
    replicate = {
        well: index
        for wells in group_wells.values()
        for index, well in enumerate(wells, 1)
    }
    plate_id = plate_label or re.sub(
        r"_?96[_ -]?well[_ -]?plate[_ -]?template$",
        "",
        paths["template"].stem,
        flags=re.IGNORECASE,
    ).rstrip("_ -")
    plate_id = plate_id or paths["template"].stem
    columns = [
        "Image_ID",
        "Alignment_Source_Row",
        "Thickness_Source_Row",
        "Plate_ID",
        "Biological_Replicate_ID",
        "Group",
        "Well",
        "Plate_Row",
        "Plate_Column",
        "Technical_Replicate",
        "Image_Number",
        "Sequence_Number",
        "Alignment_Angle_Degree",
        "Alignment_File_Name",
        "Thickness_File_Name",
        "Fibronectin_Source_Row",
        "Fibronectin_File_Name",
        FN_METRIC,
        FN_THRESHOLD_COLUMN,
        FN_LOW_FLAG,
        FN_INCLUDED_FLAG,
        FN_REASON_COLUMN,
    ]
    thickness_labels = [
        f"{field} ({THICKNESS_UNITS[field]})" for field in THICKNESS_METRICS
    ]
    occupied = set(columns + thickness_labels)
    field_map, extra_values = [], {}
    for source, source_columns, records in (
        ("Alignment", alignment_columns, alignment),
        ("Thickness", thickness_columns, thickness),
        ("Fibronectin", fn_columns, fibronectin),
    ):
        for column in source_columns:
            if (
                column == "File_Name"
                or (source == "Thickness" and column in THICKNESS_METRICS)
                or (source == "Fibronectin" and column == FN_METRIC)
            ):
                continue
            output_column = (
                f"FN_Source__{column}"
                if source == "Fibronectin"
                else f"{source}_Source__{column}"
                if column in ("Image_ID", "SourceCSV")
                else column
                if column not in occupied
                else f"{source}_Source__{column}"
            )
            while output_column in occupied:
                output_column = f"{source}_Source__{output_column}"
            columns.append(output_column)
            occupied.add(output_column)
            values = (
                [row["numbers"][column] for row in records]
                if source == "Alignment" and column == metric
                else [row["raw"][column] or None for row in records]
                if column in ("Image_ID", "SourceCSV")
                else metadata_column([row["raw"][column] for row in records])
            )
            extra_values[output_column] = {
                row["image_id"]: value for row, value in zip(records, values)
            }
            field_map.append(
                {
                    "Source": source,
                    "Original_Column": column,
                    "Output_Column": output_column,
                }
            )
            if source == "Alignment" and column == metric:
                output_metric = output_column
    columns += thickness_labels
    thickness_lookup = {row["image_id"]: row for row in thickness}
    fn_lookup = {row["image_id"]: row for row in fibronectin}
    field_map.extend(
        [
            {
                "Source": "Fibronectin",
                "Original_Column": "File_Name",
                "Output_Column": "Fibronectin_File_Name",
            },
            {
                "Source": "Fibronectin",
                "Original_Column": FN_METRIC,
                "Output_Column": FN_METRIC,
            },
        ]
    )
    merged = []
    for row in alignment:
        image_id, well = row["image_id"], row["well"]
        partner = thickness_lookup[image_id]
        fn_partner = fn_lookup[image_id]
        fn_percent = fn_partner["numbers"][FN_METRIC]
        below = fn_percent < fn_threshold
        result = dict(
            zip(
                columns[:15],
                [
                    image_id,
                    row["source_row"],
                    partner["source_row"],
                    plate_id,
                    None,
                    well_map[well],
                    well,
                    well[0],
                    int(well[1:]),
                    replicate[well],
                    row["image_number"],
                    row["sequence_number"],
                    float(angle_label),
                    row["raw"]["File_Name"],
                    partner["raw"]["File_Name"],
                ],
            )
        )
        result.update(
            {
                "Fibronectin_Source_Row": fn_partner["source_row"],
                "Fibronectin_File_Name": fn_partner["raw"]["File_Name"],
                FN_METRIC: fn_percent,
                FN_THRESHOLD_COLUMN: fn_threshold,
                FN_LOW_FLAG: below,
                FN_INCLUDED_FLAG: not below,
                FN_REASON_COLUMN: f"FN_Area_Percent < {fn_threshold:g}%"
                if below
                else None,
            }
        )
        result.update(
            {
                column: values[image_id]
                for column, values in extra_values.items()
            }
        )
        result.update(
            {
                label: partner["numbers"][field]
                for label, field in zip(thickness_labels, THICKNESS_METRICS)
            }
        )
        merged.append(result)
    ranks = {group: index for index, group in enumerate(group_order)}
    merged.sort(
        key=lambda row: (
            ranks[row["Group"]],
            row["Technical_Replicate"],
            optional_number_sort(row["Sequence_Number"]),
            optional_number_sort(row["Image_Number"]),
            row["Image_ID"],
        )
    )
    retained, excluded, group_counts, well_counts = filter_counts(
        merged, group_order, group_wells
    )
    empty_groups = [
        row["Group"] for row in group_counts if not row["Retained_Images"]
    ]
    singletons = [
        row["Group"] for row in group_counts if row["Retained_Images"] == 1
    ]
    log.event(
        "INFO",
        "FN filter",
        f"FN area < {fn_threshold:g}%: {len(excluded)} flagged; "
        f"{len(retained)} retained for filtered plots; "
        f"all {len(merged)} rows preserved.",
    )
    if empty_groups:
        log.event(
            "WARNING",
            "Empty filtered groups",
            ", ".join(empty_groups) + "; positions retained with n=0.",
        )
    if singletons:
        log.event(
            "INFO",
            "Single-image filtered groups",
            ", ".join(singletons) + "; one point, no box.",
        )
    checks = [
        (
            "Overall validation status",
            "PASS",
            "All mandatory input checks passed",
        ),
        ("Alignment image count", len(alignment), "Input records"),
        ("Thickness image count", len(thickness), "Input records"),
        ("Fibronectin image count", len(fibronectin), "Input records"),
        (
            "One-to-one matched IDs",
            len(merged),
            (
                "Complete original filenames; exact alignment source-stem "
                "mapping, independent of row order"
            ),
        ),
        ("Duplicate alignment IDs", 0, "No duplicate IDs"),
        ("Duplicate thickness IDs", 0, "No duplicate IDs"),
        ("Duplicate fibronectin IDs", 0, "No duplicate IDs"),
        ("Filename parsing coverage", "100%", "All three tables"),
        (
            "Template annotation coverage",
            "100%",
            f"{matched_count}/{len(alignment)} images",
        ),
        ("Wells without groups", 0, "None"),
        ("Groups", len(group_order), "Template row-major order"),
        (
            "Result wells",
            len(image_counts),
            "Technical-replicate wells with images",
        ),
        ("Alignment metric", metric, "Detected from the original column name"),
        (
            "Analysis angle",
            angle_label,
            "Degrees; original representation preserved",
        ),
        (
            "Numeric completeness",
            "PASS",
            "No missing or non-finite required values",
        ),
        ("Alignment range", "PASS", "0-100 inclusive"),
        ("Thickness values", "PASS", "All five metrics are non-negative"),
        ("Min <= Median <= Max", "PASS", "Every thickness record"),
        (
            "FN percentage range",
            "PASS",
            "FN_Area_Percent is between 0 and 100 inclusive",
        ),
        (
            "FN percentage / pixel counts",
            "PASS" if pixel_counts_checked else "Not supplied",
            "Reconciled where both source pixel-count columns are available",
        ),
        (
            "FN area threshold (%)",
            fn_threshold,
            (
                "Exclude only values strictly below this threshold; equality "
                "is retained"
            ),
        ),
        (
            "Images in full data",
            len(merged),
            "Every matched image is preserved",
        ),
        (
            "Images in filtered data",
            len(retained),
            "FN area is greater than or equal to the threshold",
        ),
        (
            "Excluded from filtered plots",
            len(excluded),
            "Image-level FN filter only; no entire-well exclusion",
        ),
        ("Removed from full data", 0, "All source measurements retained"),
        (
            "Empty filtered groups",
            len(empty_groups),
            ", ".join(empty_groups) or "None",
        ),
        (
            "Single-image filtered groups",
            len(singletons),
            ", ".join(singletons) or "None",
        ),
        (
            "Generated plots",
            13,
            "One FN plot, six all-image plots, and six filtered plots",
        ),
        (
            "Template wells without images",
            len(unused_wells),
            ", ".join(unused_wells) or "None",
        ),
        (
            "Statistical tests",
            "Not performed",
            "No p-values or summary-statistics tables",
        ),
        (
            "Biological replicates",
            "Not assigned",
            "Biological_Replicate_ID is blank",
        ),
        (
            "Thickness units",
            "; ".join(
                f"{field}: {unit}" for field, unit in THICKNESS_UNITS.items()
            ),
            (
                "Area uses square micrometres; thickness statistics use "
                "micrometres. Values are unchanged."
            ),
        ),
    ]
    return {
        "columns": columns,
        "rows": merged,
        "metric": output_metric,
        "angle_label": angle_label,
        "retained_rows": retained,
        "excluded_rows": excluded,
        "fn_threshold": fn_threshold,
        "group_filter_counts": group_counts,
        "well_filter_counts": well_counts,
        "group_order": group_order,
        "group_wells": group_wells,
        "well_counts": dict(image_counts),
        "well_map": well_map,
        "plate_matrix": grid,
        "plate_id": plate_id,
        "template_sheet": selected,
        "thickness_units": dict(THICKNESS_UNITS),
        "qc": [
            dict(zip(["Check", "Value", "Details"], row)) for row in checks
        ],
        "field_map": field_map,
    }
