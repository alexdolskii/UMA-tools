"""
Render selected image views and optional filtered-data statistics.
"""

from __future__ import annotations

import math
from io import BytesIO
from pathlib import Path
from typing import Any

from .plot_palette import POINT_COLOR, POINT_EDGE_COLOR, report_palette
from .plot_style import (
    FONT_SIZES,
    PNG_DPI,
    panel_letter,
    plot_font,
    rc_parameters,
    short_labels,
    wrap_label,
)
from .progress import phase
from .report_schema import (
    FN_LOW_FLAG,
    FN_METRIC,
    LOW_FN_EDGE_COLOR,
    PLOT_NAMES,
    EventLogger,
    ReportData,
    ValidationError,
)


def box_definition(values):
    """
    Calculate linear/type-7 quartiles and 1.5-IQR whiskers for the plot
    only.
    """
    import numpy as np

    values = np.asarray(values, dtype=float)
    q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75], method="linear")
    spread = q3 - q1
    low = values[values >= q1 - 1.5 * spread]
    high = values[values <= q3 + 1.5 * spread]
    return {
        "q1": float(q1),
        "med": float(median),
        "q3": float(q3),
        "whislo": float(min(q1, low.min())) if len(low) else float(q1),
        "whishi": float(max(q3, high.max())) if len(high) else float(q3),
        "fliers": [],
    }


def _plot_specs(data):
    """Plot selected metrics; keep all measurements for tables and tests."""
    return [
        ("Fibronectin", FN_METRIC, "Fibronectin coverage", "%"),
        (
            "Alignment",
            data["metric"],
            f"Fibers aligned within ±{data['angle_label']}°",
            "%",
        ),
        (
            "Area",
            "Area (µm²)",
            "Thickness analysis: measured area",
            "µm²",
        ),
    ]


def _point_positions(data, max_replicates):
    """Assign stable image positions before any FN filtering."""
    import numpy as np

    groups = data["group_order"]
    # Compute positions from all images once. Surviving images keep the
    # same X position after filtering, even when another well is empty.
    x_positions = {}
    offsets = (
        np.linspace(-0.12, 0.12, max_replicates)
        if max_replicates > 1
        else [0.0]
    )
    for group_index, group in enumerate(groups, 1):
        for rep_index, well in enumerate(data["group_wells"][group]):
            records = [row for row in data["rows"] if row["Well"] == well]
            jitter = (
                np.linspace(-0.045, 0.045, len(records))
                if len(records) > 1
                else [0.0]
            )
            for row, delta in zip(records, jitter):
                x_positions[row["Image_ID"]] = float(
                    group_index + offsets[rep_index] + delta
                )
    return x_positions


def _draw_boxes(axis, rows, groups, field, styles):
    """Draw a box only when at least two observations are present."""
    boxes, box_groups, positions = [], [], []
    for position, group in enumerate(groups, 1):
        values = [row[field] for row in rows if row["Group"] == group]
        if len(values) >= 2:
            boxes.append(box_definition(values))
            box_groups.append(group)
            positions.append(position)
    if boxes:
        artists = axis.bxp(
            boxes,
            positions=positions,
            widths=0.58,
            showfliers=False,
            patch_artist=True,
            boxprops={
                "edgecolor": "#334155",
                "linewidth": 1.4,
            },
            medianprops={
                "linewidth": 2,
                "zorder": 4,
            },
            whiskerprops={"color": "#64748B"},
            capprops={"color": "#64748B"},
        )
        for group, box, median in zip(
            box_groups, artists["boxes"], artists["medians"]
        ):
            style = styles["conditions"][group]
            box.set_facecolor(style["color"])
            median.set_color(style["median_color"])
    return boxes, box_groups


def _draw_points(
    axis, data, rows, field, filtered, styles, x_positions, groups=None
):
    """
    Draw each image in neutral gray with its well marker and FN flag.
    """
    groups = groups if groups is not None else data["group_order"]
    plotted_ids, red_ids = [], []
    for group in groups:
        for well in data["group_wells"][group]:
            records = [row for row in rows if row["Well"] == well]
            if not records:
                continue
            outlined = [
                bool(row[FN_LOW_FLAG]) and not filtered for row in records
            ]
            axis.scatter(
                [x_positions[row["Image_ID"]] for row in records],
                [row[field] for row in records],
                s=52,
                color=styles["point_color"],
                marker=styles["wells"][well],
                edgecolors=[
                    LOW_FN_EDGE_COLOR if low else POINT_EDGE_COLOR
                    for low in outlined
                ],
                linewidths=[1.9 if low else 0.8 for low in outlined],
                alpha=1,
                zorder=5,
                clip_on=False,
            )
            plotted_ids.extend(row["Image_ID"] for row in records)
            red_ids.extend(
                row["Image_ID"] for row, low in zip(records, outlined) if low
            )
    return plotted_ids, red_ids


def _verify_points(
    axis,
    rows,
    plotted_ids,
    red_ids,
    expected_red,
    name,
    view,
):
    """Reconcile rendered image IDs, point counts and red outlines."""
    import numpy as np
    from matplotlib.colors import to_rgba

    expected_ids = {row["Image_ID"] for row in rows}
    if len(plotted_ids) != len(rows) or set(plotted_ids) != expected_ids:
        raise RuntimeError(
            f"{name} {view}: image IDs or point counts do not reconcile."
        )
    actual_points = sum(
        len(collection.get_offsets()) for collection in axis.collections
    )
    actual_red = sum(
        int(np.allclose(edge[:3], to_rgba(LOW_FN_EDGE_COLOR)[:3]))
        for collection in axis.collections
        for edge in collection.get_edgecolors()
    )
    if (
        actual_points != len(rows)
        or actual_red != expected_red
        or len(red_ids) != expected_red
    ):
        raise RuntimeError(
            f"{name} {view}: rendered point or red-outline "
            "counts are incorrect."
        )
    return actual_points, actual_red


def _plot_statistics(data, field, filtered):
    """Select supplied comparisons without recalculating p-values."""
    statistics = data.get("statistics")
    if statistics is None:
        return None, "Statistical tests disabled.", []
    stats_unit = statistics["unit"]
    if not filtered:
        return (
            stats_unit,
            "No tests on this full-data view; statistics use FN-filtered "
            "data only.",
            [],
        )
    test_unit = (
        "one mean per well, with equal well weights"
        if stats_unit == "well"
        else "individual images; images within a well are dependent"
    )
    note = (
        f"Two-sided Welch tests: {test_unit}. Holm adjustment across "
        "all seven metrics and treatment-versus-control comparisons "
        "within each color block. Adjusted p: * <0.05; ** <0.01; "
        "*** <0.001; **** <0.0001; ns: ≥0.05. Not tested: insufficient or "
        "undefined test data. Technical comparisons within one plate."
    )
    if field == FN_METRIC:
        note += " FN% comparisons describe retained images only."
    comparisons = [
        dict(comparison)
        for comparison in statistics["comparisons"]
        if comparison["Metric"] == field
    ]
    return stats_unit, note, comparisons


def _comparison_layout(comparisons, groups):
    """Reuse bracket levels when comparison spans do not overlap."""
    positions = {group: index for index, group in enumerate(groups, 1)}
    levels = []
    layout = []
    for comparison in comparisons:
        missing = [
            comparison[role]
            for role in ("Control", "Treatment")
            if comparison[role] not in positions
        ]
        if missing:
            if comparison["Status"] != "Not tested":
                raise RuntimeError(
                    "A tested comparison contains a condition with no "
                    "source images."
                )
            layout.append(
                {
                    **comparison,
                    "Annotation": "Not tested",
                    "Annotation_Reason": (
                        "Condition(s) have no source images and are not "
                        "plotted: " + ", ".join(missing)
                    ),
                    "Bracket_Level": None,
                    "Bracket_Left": None,
                    "Bracket_Right": None,
                }
            )
            continue
        control = positions[comparison["Control"]]
        treatment = positions[comparison["Treatment"]]
        left, right = sorted((control, treatment))
        if left == right:
            raise RuntimeError("A statistical comparison needs two groups.")
        for level, occupied in enumerate(levels):
            if all(right < start or left > end for start, end in occupied):
                occupied.append((left, right))
                break
        else:
            level = len(levels)
            levels.append([(left, right)])
        label = (
            comparison["Significance"]
            if comparison["Status"] == "Tested"
            else "Not tested"
        )
        if label not in {"*", "**", "***", "****", "ns", "Not tested"}:
            raise RuntimeError("Unexpected statistical plot annotation.")
        layout.append(
            {
                **comparison,
                "Annotation": label,
                "Bracket_Level": level,
                "Bracket_Left": left,
                "Bracket_Right": right,
            }
        )
    return layout, len(levels)


def _draw_comparisons(axis, comparisons, level_count):
    """Keep comparison brackets separate from calibrated data axes."""
    axis.set_axis_off()
    axis.set_ylim(0, level_count + 0.25)
    for comparison in comparisons:
        if comparison["Bracket_Level"] is None:
            continue
        left = comparison["Bracket_Left"]
        right = comparison["Bracket_Right"]
        bottom = comparison["Bracket_Level"] + 0.12
        top = bottom + 0.23
        axis.plot(
            [left, left, right, right],
            [bottom, top, top, bottom],
            color="#334155",
            linewidth=1,
        )
        axis.text(
            (left + right) / 2,
            top + 0.04,
            comparison["Annotation"],
            ha="center",
            va="bottom",
            fontsize=FONT_SIZES["axis"],
            fontweight="bold",
            color="#334155",
        )


def prepare_plot_design(data, template, log):
    """Read descriptive colors without enabling statistical tests."""
    statistics = data.get("statistics")
    if statistics is not None:
        data["plot_design"] = statistics["design"]
        data["template_theme"] = statistics.get("template_theme")
        data["template_palette"] = statistics.get("template_palette")
        return
    import openpyxl

    from .report_statistics import _color_fields, _reject_conditional_styles

    workbook = openpyxl.load_workbook(
        template, data_only=False, rich_text=True
    )
    from openpyxl.cell.rich_text import CellRichText

    design = []
    try:
        sheet = workbook[data["template_sheet"]]
        data["template_theme"] = workbook.loaded_theme
        data["template_palette"] = list(workbook._colors)
        try:
            _reject_conditional_styles(sheet)
        except ValidationError:
            log.event(
                "INFO",
                "Plot panels",
                "Conditional plate styles cannot define panels; "
                "using all conditions together. Statistics remain disabled.",
            )
            data["plot_design"] = []
            return
        for row, letter in enumerate("ABCDEFGH", 2):
            for column in range(1, 13):
                well = f"{letter}{column:02d}"
                if well not in data["well_map"]:
                    continue
                cell = sheet.cell(row, column + 1)
                try:
                    fields = _color_fields(cell)
                except ValueError:
                    fields = {"Color_Code": None}
                design.append(
                    {
                        "Group": data["well_map"][well],
                        "Well": well,
                        "Excel_Cell": cell.coordinate,
                        "Is_Control": None
                        if isinstance(cell.value, CellRichText)
                        else bool(cell.font.bold),
                        **fields,
                    }
                )
    finally:
        workbook.close()
    data["plot_design"] = design


def plot_panels(data):
    """Partition conditions once, retaining empty filtered groups."""
    groups = data["group_order"]
    statistics = data.get("statistics") or {}
    design = data.get("plot_design", statistics.get("design", []))
    by_well = {row["Well"]: row for row in design}
    group_colors = {}
    for group in groups:
        wells = (
            [well for well, name in data["well_map"].items() if name == group]
            if data.get("well_map")
            else data["group_wells"][group]
        )
        colors = {by_well.get(well, {}).get("Color_Code") for well in wells}
        # Ambiguous or incomplete descriptive markup must not split a
        # condition, discard observations, or require controls.
        if len(colors) != 1 or None in colors:
            group_colors = {}
            break
        group_colors[group] = colors.pop()
    blocks = {}
    for group in groups:
        blocks.setdefault(group_colors.get(group, ""), []).append(group)
    if any(
        row["Control"] in group_colors
        and row["Treatment"] in group_colors
        and group_colors[row["Control"]] != group_colors[row["Treatment"]]
        for row in statistics.get("comparisons", [])
    ):
        blocks = {"": groups}
    context, _ = short_labels(groups)
    panels = []
    for index, (color, members) in enumerate(blocks.items(), 1):
        prefix, labels = short_labels(members)
        heading = prefix
        if context and prefix.startswith(context):
            heading = prefix[len(context) :].strip()
        panels.append(
            {
                "id": panel_letter(index),
                "color_code": color,
                "groups": members,
                "labels": labels,
                "context": context,
                "title": heading
                or (f"Block {index}" if color else "All conditions"),
            }
        )
    return panels


def _legend_handles(name, threshold, filtered, markers):
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            linestyle="none",
            markerfacecolor=POINT_COLOR,
            markeredgecolor=POINT_EDGE_COLOR,
            markersize=8,
            label=f"Technical well {index}",
        )
        for index, marker in enumerate(markers, 1)
    ]
    if not filtered:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=POINT_COLOR,
                markeredgecolor=LOW_FN_EDGE_COLOR,
                markeredgewidth=1.9,
                markersize=8,
                label=f"Red outline: FN < {threshold:g}%",
            )
        )
    if name == "Fibronectin":
        handles.append(
            Line2D(
                [0],
                [0],
                color=LOW_FN_EDGE_COLOR,
                linestyle="--",
                label=f"FN cutoff: {threshold:g}%",
            )
        )
    return handles


def _render_plot(
    data,
    directory,
    spec,
    filtered,
    styles,
    shared_upper,
    x_positions,
    log,
    plot_format="pdf",
):
    """Render colored boxes and gray image points with well shapes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    name, field, base_title, unit = spec
    groups, threshold = data["group_order"], data["fn_threshold"]
    rows = data["retained_rows"] if filtered else data["rows"]
    panels = plot_panels(data)
    stats_unit, stats_note, supplied = _plot_statistics(data, field, filtered)
    # Keep untestable comparisons involving unimaged conditions in the
    # manifest and Statistics table; there is no invented X position.
    all_comparisons, _ = _comparison_layout(supplied, groups)
    missing = [r for r in all_comparisons if r["Bracket_Level"] is None]
    if missing:
        stats_note += (
            " Comparisons involving conditions with no source images "
            "are listed in the Statistics sheet."
        )
    counts = {g: sum(r["Group"] == g for r in rows) for g in groups}
    well_counts = {
        g: len({r["Well"] for r in rows if r["Group"] == g}) for g in groups
    }
    columns = min(2, len(panels))
    panel_width = max(7.2, 1.85 * max(len(p["groups"]) for p in panels))
    width = columns * panel_width
    view = f"FN coverage ≥ {threshold:g}%" if filtered else "All images"
    title = f"{base_title} — {view}"
    heading = wrap_label(
        title, (width - 0.6) * 72, FONT_SIZES["title"], "bold"
    )
    subtitle = wrap_label(
        panels[0]["context"], (width - 0.6) * 72, FONT_SIZES["panel"]
    )
    handles = _legend_handles(name, threshold, filtered, styles["markers"])
    legend_columns = min(4, max(2, int(width / 3)))
    legend_rows = math.ceil(len(handles) / legend_columns)
    panel_labels, layouts, local_positions = {}, {}, {}
    for panel in panels:
        label_width = (
            (panel_width - 1.3) * 72 / (len(panel["groups"]) + 0.2) * 0.9
        )
        panel_labels[panel["id"]] = {
            g: wrap_label(label, label_width, FONT_SIZES["axis"])
            for g, label in panel["labels"].items()
        }
        comparisons = [
            r
            for r in supplied
            if r["Control"] in panel["groups"]
            and r["Treatment"] in panel["groups"]
        ]
        layouts[panel["id"]] = _comparison_layout(comparisons, panel["groups"])
        for local_index, group in enumerate(panel["groups"], 1):
            shift = local_index - (groups.index(group) + 1)
            for row in data["rows"]:
                if row["Group"] == group:
                    image_id = row["Image_ID"]
                    local_positions[image_id] = x_positions[image_id] + shift
    label_lines = max(
        s.count("\n") + 1
        for labels in panel_labels.values()
        for s in labels.values()
    )
    levels = max(levels for _, levels in layouts.values())
    mask_caption = data.get("fn_mask_settings", {}).get(
        "caption", "FN mask intensity thresholds: not recorded."
    )
    caption = f"{mask_caption}\nCoverage filter: FN ≥ {threshold:g}%; "
    caption += (
        f"{len(rows)} retained images. "
        if filtered
        else f"all {len(rows)} images shown; lower coverage outlined in red. "
    )
    caption += (
        "\nOne point = one image. Box: median and IQR; whiskers: 1.5 IQR. "
    )
    caption += (
        f"Welch + Holm; test unit: {stats_unit}. Symbols: Statistics table."
        if filtered and stats_unit
        else "Descriptive view; no tests."
    )
    caption += (
        "\nBox color = condition; gray points = images; "
        "shape = technical well within condition."
    )
    if any(style["is_control"] for style in styles["conditions"].values()):
        caption += " Lavender = bold-marked control."
    display_caption = wrap_label(
        caption, (width - 0.6) * 72, FONT_SIZES["note"]
    )
    source = data.get("source_name", data.get("plate_id", ""))
    footer_label = wrap_label(
        f"{source} · Run: {data.get('run_id', '')}".strip(" ·"),
        (width - 0.6) * 72,
        FONT_SIZES["note"],
    )
    header_height = (
        0.4 * (heading.count("\n") + 1)
        + (0.32 * (subtitle.count("\n") + 1) if subtitle else 0)
        + 0.35 * legend_rows
        + 0.2
    )
    body_height = math.ceil(len(panels) / columns) * (
        4.2 + 0.28 * label_lines + 0.45 * levels + 0.55
    )
    footer_height = (
        0.23 * (display_caption.count("\n") + footer_label.count("\n") + 2)
        + 0.32
    )
    height = header_height + body_height + footer_height + 0.4
    stem = (
        "fibronectin_boxplot"
        if name == "Fibronectin"
        else "alignment_boxplot"
        if name == "Alignment"
        else f"thickness_{name.lower()}_boxplot"
    ) + ("_filtered" if filtered else "")
    plot_id = name + (" Filtered" if filtered else " Plot")
    png_path, pdf_path = (
        directory / (stem + ".png"),
        directory / (stem + ".pdf"),
    )
    rendered_ids, red_ids, boxes_by_group, annotations = [], [], {}, []
    with plt.rc_context(rc_parameters()):
        fig = plt.figure(figsize=(width, height), layout="constrained")
        fig.get_layout_engine().set(
            w_pad=0.12, h_pad=0.12, hspace=0.06, wspace=0.04
        )
        try:
            outer = fig.add_gridspec(
                3, 1, height_ratios=[header_height, body_height, footer_height]
            )
            header = fig.add_subplot(outer[0])
            header.set_axis_off()
            header.set_label("header")
            header.text(
                0.5,
                1,
                heading,
                ha="center",
                va="top",
                fontsize=FONT_SIZES["title"],
                fontweight="bold",
                transform=header.transAxes,
            )
            if subtitle:
                header.text(
                    0.5,
                    (0.35 * legend_rows + 0.12) / header_height,
                    subtitle,
                    ha="center",
                    va="bottom",
                    fontsize=FONT_SIZES["panel"],
                    transform=header.transAxes,
                )
            header.legend(
                handles=handles,
                loc="lower center",
                ncol=legend_columns,
                frameon=False,
                fontsize=FONT_SIZES["note"],
                borderaxespad=0,
            )
            grid = outer[1].subgridspec(
                math.ceil(len(panels) / columns), columns
            )
            for index, panel in enumerate(panels):
                local_rows = [r for r in rows if r["Group"] in panel["groups"]]
                local_layout, _ = layouts[panel["id"]]
                slot = grid[index // columns, index % columns]
                if levels:
                    subgrid = slot.subgridspec(
                        2, 1, height_ratios=[0.45 * levels, 4.2]
                    )
                    bracket_axis = fig.add_subplot(subgrid[0])
                    bracket_axis.set_label("comparisons_" + panel["id"])
                    axis = fig.add_subplot(subgrid[1], sharex=bracket_axis)
                    _draw_comparisons(bracket_axis, local_layout, levels)
                    title_axis = bracket_axis
                else:
                    axis = fig.add_subplot(slot)
                    title_axis = axis
                axis.set_label("data_" + panel["id"])
                panel_title = wrap_label(
                    f"{panel['id']}  {panel['title']}",
                    (panel_width - 1.3) * 72,
                    FONT_SIZES["panel"],
                    "bold",
                )
                title_axis.set_title(
                    panel_title,
                    fontsize=FONT_SIZES["panel"],
                    fontweight="bold",
                    pad=12,
                )
                boxes, box_groups = _draw_boxes(
                    axis, local_rows, panel["groups"], field, styles
                )
                ids, reds = _draw_points(
                    axis,
                    data,
                    local_rows,
                    field,
                    filtered,
                    styles,
                    local_positions,
                    panel["groups"],
                )
                _verify_points(
                    axis,
                    local_rows,
                    ids,
                    reds,
                    0 if filtered else sum(r[FN_LOW_FLAG] for r in local_rows),
                    name,
                    view,
                )
                rendered_ids.extend(ids)
                red_ids.extend(reds)
                boxes_by_group.update(zip(box_groups, boxes))
                annotations.extend(
                    {**r, "Panel": panel["id"]} for r in local_layout
                )
                axis.set_ylim(0, shared_upper[field])
                axis.set_xlim(0.4, len(panel["groups"]) + 0.6)
                ylabel = (
                    "Fibronectin coverage (%)"
                    if name == "Fibronectin"
                    else "Fibers aligned (%)"
                    if name == "Alignment"
                    else f"{base_title} ({unit})"
                )
                axis.set_ylabel(ylabel, fontsize=FONT_SIZES["axis"])
                labels = panel_labels[panel["id"]]
                axis.set_xticks(
                    range(1, len(panel["groups"]) + 1),
                    [labels[g] for g in panel["groups"]],
                    fontsize=FONT_SIZES["axis"],
                )
                lines = max(s.count("\n") + 1 for s in labels.values())
                for position, group in enumerate(panel["groups"], 1):
                    image_word = "image" if counts[group] == 1 else "images"
                    well_word = "well" if well_counts[group] == 1 else "wells"
                    sample = (
                        f"n={counts[group]} {image_word}\n"
                        f"{well_counts[group]} {well_word}"
                    )
                    label = axis.annotate(
                        sample,
                        (position, 0),
                        xycoords=axis.get_xaxis_transform(),
                        xytext=(0, -(12 + lines * FONT_SIZES["axis"] * 1.2)),
                        textcoords="offset points",
                        ha="center",
                        va="top",
                        fontsize=FONT_SIZES["sample"],
                        annotation_clip=False,
                    )
                    label.set_gid("sample_size")
                axis.tick_params(axis="x", pad=6)
                axis.tick_params(axis="y", labelsize=FONT_SIZES["axis"])
                axis.ticklabel_format(axis="y", style="sci", scilimits=(-4, 5))
                axis.yaxis.get_offset_text().set_fontsize(FONT_SIZES["sample"])
                axis.grid(axis="y", color="#D7DEE8", linewidth=0.8)
                axis.set_axisbelow(True)
                axis.spines[["top", "right"]].set_visible(False)
                axis.spines[["left", "bottom"]].set_color("#94A3B8")
                if name == "Fibronectin":
                    axis.axhline(
                        threshold,
                        color=LOW_FN_EDGE_COLOR,
                        linestyle="--",
                        linewidth=1.3,
                    )
                if not local_rows:
                    axis.text(
                        0.5,
                        0.5,
                        "No retained images",
                        ha="center",
                        va="center",
                        transform=axis.transAxes,
                        fontsize=FONT_SIZES["note"],
                    )
            if len(rendered_ids) != len(rows) or set(rendered_ids) != {
                r["Image_ID"] for r in rows
            }:
                raise RuntimeError(
                    "Panel assignment changed the image population"
                )
            footer = fig.add_subplot(outer[2])
            footer.set_axis_off()
            footer.set_label("footer")
            footer.text(
                0,
                1,
                display_caption,
                fontsize=FONT_SIZES["note"],
                va="top",
                transform=footer.transAxes,
            )
            footer.text(
                0,
                0,
                footer_label,
                fontsize=FONT_SIZES["note"],
                va="bottom",
                transform=footer.transAxes,
            )
            if plot_format in ("png", "both"):
                fig.savefig(png_path, dpi=PNG_DPI, facecolor="white")
            if plot_format in ("pdf", "both"):
                fig.savefig(pdf_path, facecolor="white")
            if plot_format == "pdf":
                with BytesIO() as preview:
                    fig.savefig(
                        preview, format="png", dpi=PNG_DPI, facecolor="white"
                    )
                    data.setdefault("plot_previews", {})[plot_id] = (
                        preview.getvalue()
                    )
        finally:
            plt.close(fig)
    # Stable order in metadata, even if blocks change panel placement.
    order = {
        r["Image_ID"]: (
            groups.index(r["Group"]),
            data["group_wells"][r["Group"]].index(r["Well"]),
            i,
        )
        for i, r in enumerate(data["rows"])
    }
    rendered_ids.sort(key=order.get)
    red_ids.sort(key=order.get)
    plot = {
        "name": name,
        "plot_id": plot_id,
        "sheet": "Plots",
        "view": "Filtered" if filtered else "All images",
        "metric": field,
        "title": title,
        "unit": unit,
        "path": str(pdf_path if plot_format == "pdf" else png_path),
        "png_file": str(png_path) if plot_format in ("png", "both") else "",
        "pdf_file": str(pdf_path) if plot_format in ("pdf", "both") else "",
        "format": plot_format,
        "font": plot_font(),
        "png_dpi": PNG_DPI,
        "caption": caption,
        "fn_mask_caption": mask_caption,
        "statistics_note": stats_note,
        "point_count": len(rendered_ids),
        "red_outline_count": len(red_ids),
        "y_min": 0.0,
        "y_max": shared_upper[field],
        "width": width,
        "height": height,
        "panels": panels,
        "group_order": groups,
        "group_counts": counts,
        "group_well_counts": well_counts,
        "statistics_unit": stats_unit,
        "statistical_comparisons": annotations + missing,
        "box_groups": [g for g in groups if g in boxes_by_group],
        "empty_groups": [g for g in groups if counts[g] == 0],
        "singleton_groups": [g for g in groups if counts[g] == 1],
        "plotted_image_ids": rendered_ids,
        "red_outline_image_ids": red_ids,
        "x_positions": {key: x_positions[key] for key in rendered_ids},
        "rendered_x_positions": {
            key: local_positions[key] for key in rendered_ids
        },
        "condition_styles": styles["conditions"],
        "point_color": styles["point_color"],
        "point_edge_color": styles["point_edge_color"],
        "well_markers": styles["wells"],
        "technical_well_markers": styles["markers"],
        "palette_blocks": styles["blocks"],
        "boxes": boxes_by_group,
    }
    log.event(
        "PASS",
        "Plot",
        f"{plot_id}: {len(panels)} panel(s); "
        f"{len(rendered_ids)} points; {len(red_ids)} red outlines; "
        f"Y=0 to {shared_upper[field]:.6g} {unit}",
    )
    return plot


def create_plots(
    data: ReportData,
    directory: Path,
    log: EventLogger,
    plot_format="pdf",
) -> list[dict[str, Any]]:
    """Build full and filtered views, independently of export format."""
    if plot_format not in ("pdf", "png", "both"):
        raise ValueError("plot_format must be pdf, png, or both")
    directory.mkdir()
    max_replicates = max(map(len, data["group_wells"].values()))
    styles = report_palette(data, plot_panels(data), log)
    data["plot_palette"] = styles
    specs = _plot_specs(data)
    shared_upper = {
        field: (
            100.0
            if name in ("Alignment", "Fibronectin")
            else max(row[field] for row in data["rows"]) * 1.12 or 1.0
        )
        for name, field, _, _ in specs
    }
    x_positions = _point_positions(data, max_replicates)
    jobs = [(spec, False) for spec in specs] + [(spec, True) for spec in specs]
    plots = []
    for index, (spec, filtered) in enumerate(jobs, 1):
        phase(f"Plot {index}/{len(jobs)}: {spec[0]}, filtered={filtered}")
        plots.append(
            _render_plot(
                data,
                directory,
                spec,
                filtered,
                styles,
                shared_upper,
                x_positions,
                log,
                plot_format,
            )
        )
    if [plot["plot_id"] for plot in plots] != PLOT_NAMES:
        raise RuntimeError("The required plot order was not preserved.")
    return plots
