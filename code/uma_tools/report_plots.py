"""
Render fourteen image views and optional filtered-data statistics.
"""

from __future__ import annotations

import colorsys
import math
import textwrap
from pathlib import Path
from typing import Any

from .report_schema import (
    BASE_COLORS,
    FN_LOW_FLAG,
    FN_METRIC,
    LOW_FN_EDGE_COLOR,
    SHEET_NAMES,
    THICKNESS_METRICS,
    THICKNESS_UNITS,
    EventLogger,
    ReportData,
)


def replicate_colors(count):
    colors = BASE_COLORS.copy()
    index = 0
    while len(colors) < count:
        rgb = colorsys.hls_to_rgb((index * 0.61803398875) % 1, 0.44, 0.62)
        color = "#" + "".join(f"{round(value * 255):02x}" for value in rgb)
        if color.lower() not in [item.lower() for item in colors]:
            colors.append(color)
        index += 1
    return colors[:count]


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
    """Keep the original full-data metric order and unit labels."""
    specs = [
        ("Fibronectin", FN_METRIC, "Fibronectin-positive area by group", "%"),
        (
            "Alignment",
            data["metric"],
            f"Fibers aligned within ±{data['angle_label']}° by group",
            "%",
        ),
    ]
    specs += [
        (
            field,
            f"{field} ({THICKNESS_UNITS[field]})",
            f"{field} by group",
            THICKNESS_UNITS[field],
        )
        for field in THICKNESS_METRICS
    ]
    return specs


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


def _draw_boxes(axis, rows, groups, field):
    """Draw a box only when at least two observations are present."""
    boxes, box_groups, positions = [], [], []
    for position, group in enumerate(groups, 1):
        values = [row[field] for row in rows if row["Group"] == group]
        if len(values) >= 2:
            boxes.append(box_definition(values))
            box_groups.append(group)
            positions.append(position)
    if boxes:
        axis.bxp(
            boxes,
            positions=positions,
            widths=0.58,
            showfliers=False,
            patch_artist=True,
            boxprops={
                "facecolor": "#E8EDF3",
                "edgecolor": "#334155",
                "linewidth": 1.4,
            },
            medianprops={
                "color": "#111827",
                "linewidth": 2,
                "zorder": 4,
            },
            whiskerprops={"color": "#64748B"},
            capprops={"color": "#64748B"},
        )
    return boxes, box_groups


def _draw_points(axis, data, rows, field, filtered, colors, x_positions):
    """
    Draw each image once using its original well color and FN flag.
    """
    groups = data["group_order"]
    plotted_ids, red_ids = [], []
    for group in groups:
        for rep_index, well in enumerate(data["group_wells"][group]):
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
                color=colors[rep_index],
                edgecolors=[
                    LOW_FN_EDGE_COLOR if low else "white" for low in outlined
                ],
                linewidths=[1.9 if low else 0.6 for low in outlined],
                alpha=0.95,
                zorder=3,
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
        "*** <0.001; ns: ≥0.05. Not tested: insufficient or "
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
        if label not in {"*", "**", "***", "ns", "Not tested"}:
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
            fontsize=9,
            color="#334155",
        )


def _render_plot(
    data,
    directory,
    spec,
    filtered,
    colors,
    shared_upper,
    x_positions,
    log,
):
    """Render, verify and save one complete plot, always closing it."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    name, field, base_title, unit = spec
    threshold, groups = data["fn_threshold"], data["group_order"]
    max_replicates = len(colors)
    rows = data["retained_rows"] if filtered else data["rows"]
    stats_unit, stats_note, comparisons = _plot_statistics(
        data, field, filtered
    )
    comparisons, comparison_levels = _comparison_layout(comparisons, groups)
    if any(row["Bracket_Level"] is None for row in comparisons):
        stats_note += (
            " Comparisons involving conditions with no source images "
            "are listed in the Statistics sheet."
        )
    counts = {
        group: sum(row["Group"] == group for row in rows) for group in groups
    }
    well_counts = {
        group: len({row["Well"] for row in rows if row["Group"] == group})
        for group in groups
    }
    labels = [
        textwrap.fill(
            group, width=26, break_long_words=True, break_on_hyphens=False
        )
        + (
            f"\nn_images={counts[group]}; n_wells={well_counts[group]}"
            if stats_unit is not None
            else f"\nn={counts[group]}"
        )
        for group in groups
    ]
    width = max(14.2, 1.2 * len(labels))
    legend_count = (
        max_replicates
        + (0 if filtered else 1)
        + (1 if name == "Fibronectin" else 0)
    )
    legend_rows = math.ceil(legend_count / 4)
    height = (
        8.7
        + 0.28 * (legend_rows - 1)
        + 0.2 * max(label.count("\n") for label in labels)
    )
    annotation_height = 0.38 * comparison_levels
    caption_height = 0.58 if stats_unit is not None else 0
    base_height = height
    height += annotation_height + caption_height
    view = f"FN area ≥ {threshold:g}%" if filtered else "All images"
    title = f"{base_title} — {view}"
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "text.usetex": False,
            "text.parse_math": False,
        }
    ):
        if comparison_levels:
            fig, (comparison_axis, axis) = plt.subplots(
                2,
                1,
                figsize=(width, height),
                dpi=160,
                sharex=True,
                gridspec_kw={
                    "height_ratios": [annotation_height, base_height],
                },
            )
        else:
            fig, axis = plt.subplots(figsize=(width, height), dpi=160)
        try:
            boxes, box_groups = _draw_boxes(axis, rows, groups, field)
            plotted_ids, red_ids = _draw_points(
                axis,
                data,
                rows,
                field,
                filtered,
                colors,
                x_positions,
            )
            actual_points, actual_red = _verify_points(
                axis,
                rows,
                plotted_ids,
                red_ids,
                0 if filtered else len(data["excluded_rows"]),
                name,
                view,
            )
            upper = shared_upper[field]
            axis.set_ylim(0, upper)
            axis.set_xlim(0.4, len(labels) + 0.6)
            if comparison_levels:
                _draw_comparisons(
                    comparison_axis, comparisons, comparison_levels
                )
            axis.set_ylabel(
                "Fibronectin-positive area (%)"
                if name == "Fibronectin"
                else "Fibers aligned (%)"
                if name == "Alignment"
                else field
            )
            axis.set_xlabel("Group")
            axis.set_xticks(range(1, len(labels) + 1))
            axis.set_xticklabels(labels, rotation=32, ha="right", fontsize=10)
            axis.grid(axis="y", color="#D7DEE8", linewidth=0.8)
            axis.set_axisbelow(True)
            axis.spines[["top", "right"]].set_visible(False)
            axis.spines[["left", "bottom"]].set_color("#94A3B8")
            if not rows:
                axis.text(
                    0.5,
                    0.5,
                    f"No images meet FN area ≥ {threshold:g}%",
                    transform=axis.transAxes,
                    ha="center",
                    va="center",
                    color="#475569",
                )
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="none",
                    markerfacecolor=colors[index],
                    markeredgecolor="white",
                    markersize=8,
                    label=f"Technical replicate {index + 1}",
                )
                for index in range(max_replicates)
            ]
            if not filtered:
                handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        linestyle="none",
                        markerfacecolor="#CBD5E1",
                        markeredgecolor=LOW_FN_EDGE_COLOR,
                        markeredgewidth=1.9,
                        markersize=8,
                        label=f"Red outline: FN area < {threshold:g}%",
                    )
                )
            if name == "Fibronectin":
                axis.axhline(
                    threshold,
                    color=LOW_FN_EDGE_COLOR,
                    linestyle="--",
                    linewidth=1.3,
                )
                handles.append(
                    Line2D(
                        [0],
                        [0],
                        color=LOW_FN_EDGE_COLOR,
                        linestyle="--",
                        label=f"FN area threshold: {threshold:g}%",
                    )
                )
            fig.suptitle(title, fontsize=16, y=0.975)
            fig.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.93),
                ncol=min(4, legend_count),
                frameon=False,
                fontsize=10,
            )
            count_note = (
                (
                    f"{len(rows)} retained; "
                    f"{len(data['excluded_rows'])} excluded "
                    "from this view."
                )
                if filtered
                else (
                    f"{len(rows)} images; "
                    f"{len(data['excluded_rows'])} below the FN threshold."
                )
            )
            caption = (
                count_note
                + " Point fill identifies the original technical-replicate "
                "well. Points and boxes always represent images.\n"
                "Boxes: 1.5×IQR whiskers; n=1: point only; n=0: no point "
                "or box.\n"
                + textwrap.fill(stats_note, width=max(100, int(width * 13)))
            )
            fig.text(
                0.5,
                0.026,
                caption,
                ha="center",
                fontsize=9,
                color="#475569",
            )
            fig.tight_layout(
                h_pad=0.2,
                rect=(
                    0.015,
                    0.08 + caption_height / height,
                    0.985,
                    0.90 - 0.03 * (legend_rows - 1),
                ),
            )
            stem = (
                "fibronectin_boxplot"
                if name == "Fibronectin"
                else "alignment_boxplot"
                if name == "Alignment"
                else f"thickness_{name.lower()}_boxplot"
            )
            path = directory / (
                stem + ("_filtered" if filtered else "") + ".png"
            )
            fig.savefig(path, facecolor="white")
            plot = {
                "name": name,
                "sheet": name + (" Filtered" if filtered else " Plot"),
                "view": "Filtered" if filtered else "All images",
                "metric": field,
                "title": title,
                "unit": unit,
                "path": str(path),
                "point_count": actual_points,
                "red_outline_count": actual_red,
                "y_min": 0.0,
                "y_max": upper,
                "width": width,
                "height": height,
                "group_order": groups,
                "group_counts": counts,
                "group_well_counts": well_counts,
                "statistics_unit": stats_unit,
                "statistics_note": stats_note,
                "statistical_comparisons": comparisons,
                "box_groups": box_groups,
                "empty_groups": [
                    group for group in groups if counts[group] == 0
                ],
                "singleton_groups": [
                    group for group in groups if counts[group] == 1
                ],
                "plotted_image_ids": plotted_ids,
                "red_outline_image_ids": red_ids,
                "x_positions": {
                    image_id: x_positions[image_id] for image_id in plotted_ids
                },
                "technical_replicate_colors": colors,
                "boxes": {group: box for group, box in zip(box_groups, boxes)},
            }
            log.event(
                "PASS",
                "Plot",
                f"{plot['sheet']}: {actual_points} points; "
                f"{actual_red} red outlines; "
                f"Y=0 to {upper:.6g} {unit}",
            )
        finally:
            plt.close(fig)
    return plot


def create_plots(
    data: ReportData, directory: Path, log: EventLogger
) -> list[dict[str, Any]]:
    """
    Build seven full-data and seven filtered plots with stable colors
    and axes.
    """
    directory.mkdir()
    max_replicates = max(map(len, data["group_wells"].values()))
    colors = replicate_colors(max_replicates)
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
    for (name, field, base_title, unit), filtered in jobs:
        plots.append(
            _render_plot(
                data,
                directory,
                (name, field, base_title, unit),
                filtered,
                colors,
                shared_upper,
                x_positions,
                log,
            )
        )
    if [plot["sheet"] for plot in plots] != SHEET_NAMES[:14]:
        raise RuntimeError("The required 14-plot order was not preserved.")
    return plots
