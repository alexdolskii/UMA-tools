"""
Render the original 13 plot views with stable points, colors, and axes.
"""

from __future__ import annotations

import colorsys
import math
import textwrap
from pathlib import Path
from typing import Any

from .constants import (
    BASE_COLORS,
    FN_LOW_FLAG,
    FN_METRIC,
    LOW_FN_EDGE_COLOR,
    SHEET_NAMES,
    THICKNESS_METRICS,
    THICKNESS_UNITS,
)
from .models import EventLogger, ReportData


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


def create_plots(
    data: ReportData, directory: Path, log: EventLogger
) -> list[dict[str, Any]]:
    """
    Build seven full-data and six filtered plots with stable colors and
    axes.
    """
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D

    directory.mkdir()
    threshold = data["fn_threshold"]
    groups = data["group_order"]
    max_replicates = max(map(len, data["group_wells"].values()))
    colors = replicate_colors(max_replicates)
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
    shared_upper = {
        field: (
            100.0
            if name in ("Alignment", "Fibronectin")
            else max(row[field] for row in data["rows"]) * 1.12 or 1.0
        )
        for name, field, _, _ in specs
    }
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
    jobs = [(spec, False) for spec in specs] + [
        (spec, True) for spec in specs[1:]
    ]
    plots = []
    for (name, field, base_title, unit), filtered in jobs:
        rows = data["retained_rows"] if filtered else data["rows"]
        counts = {
            group: sum(row["Group"] == group for row in rows)
            for group in groups
        }
        labels = [
            textwrap.fill(
                group, width=26, break_long_words=True, break_on_hyphens=False
            )
            + f"\nn={counts[group]}"
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
            fig, axis = plt.subplots(figsize=(width, height), dpi=160)
            try:
                boxes, box_groups, positions = [], [], []
                for position, group in enumerate(groups, 1):
                    values = [
                        row[field] for row in rows if row["Group"] == group
                    ]
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
                plotted_ids, red_ids = [], []
                for group in groups:
                    for rep_index, well in enumerate(
                        data["group_wells"][group]
                    ):
                        records = [row for row in rows if row["Well"] == well]
                        if not records:
                            continue
                        outlined = [
                            bool(row[FN_LOW_FLAG]) and not filtered
                            for row in records
                        ]
                        axis.scatter(
                            [x_positions[row["Image_ID"]] for row in records],
                            [row[field] for row in records],
                            s=52,
                            color=colors[rep_index],
                            edgecolors=[
                                LOW_FN_EDGE_COLOR if low else "white"
                                for low in outlined
                            ],
                            linewidths=[
                                1.9 if low else 0.6 for low in outlined
                            ],
                            alpha=0.95,
                            zorder=3,
                            clip_on=False,
                        )
                        plotted_ids.extend(row["Image_ID"] for row in records)
                        red_ids.extend(
                            row["Image_ID"]
                            for row, low in zip(records, outlined)
                            if low
                        )
                expected_ids = {row["Image_ID"] for row in rows}
                if (
                    len(plotted_ids) != len(rows)
                    or set(plotted_ids) != expected_ids
                ):
                    raise RuntimeError(
                        f"{name} {view}: image IDs or point counts "
                        "do not reconcile."
                    )
                actual_points = sum(
                    len(collection.get_offsets())
                    for collection in axis.collections
                )
                actual_red = sum(
                    int(np.allclose(edge[:3], to_rgba(LOW_FN_EDGE_COLOR)[:3]))
                    for collection in axis.collections
                    for edge in collection.get_edgecolors()
                )
                expected_red = 0 if filtered else len(data["excluded_rows"])
                if (
                    actual_points != len(rows)
                    or actual_red != expected_red
                    or len(red_ids) != expected_red
                ):
                    raise RuntimeError(
                        f"{name} {view}: rendered point or red-outline "
                        "counts are incorrect."
                    )
                upper = shared_upper[field]
                axis.set_ylim(0, upper)
                axis.set_xlim(0.4, len(labels) + 0.6)
                axis.set_ylabel(
                    "Fibronectin-positive area (%)"
                    if name == "Fibronectin"
                    else "Fibers aligned (%)"
                    if name == "Alignment"
                    else field
                )
                axis.set_xlabel("Group")
                axis.set_xticks(range(1, len(labels) + 1))
                axis.set_xticklabels(
                    labels, rotation=32, ha="right", fontsize=10
                )
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
                fig.text(
                    0.5,
                    0.026,
                    count_note
                    + (
                        " Point fill identifies the original "
                        "technical-replicate well.\n"
                        "Boxes: 1.5×IQR whiskers; n=1: point only; n=0: no "
                        "point or box. No statistical tests."
                    ),
                    ha="center",
                    fontsize=9,
                    color="#475569",
                )
                fig.tight_layout(
                    rect=(0.015, 0.08, 0.985, 0.90 - 0.03 * (legend_rows - 1))
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
                        image_id: x_positions[image_id]
                        for image_id in plotted_ids
                    },
                    "technical_replicate_colors": colors,
                    "boxes": {
                        group: box for group, box in zip(box_groups, boxes)
                    },
                }
                plots.append(plot)
                log.event(
                    "PASS",
                    "Plot",
                    f"{plot['sheet']}: {actual_points} points; "
                    f"{actual_red} red outlines; "
                    f"Y=0 to {upper:.6g} {unit}",
                )
            finally:
                plt.close(fig)
    if [plot["sheet"] for plot in plots] != SHEET_NAMES[:13]:
        raise RuntimeError("The required 13-plot order was not preserved.")
    return plots
