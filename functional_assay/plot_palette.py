"""Shared report styles assigned from the full plate, before exclusions."""

from __future__ import annotations

import math

from uma_tools.files import save_json
from uma_tools.plot_palette import (
    POINT_COLOR,
    POINT_EDGE_COLOR,
    condition_palette,
    well_markers,
)

PALETTE_COLUMNS = (
    "Comparison_Block",
    "Group",
    "Is_Control",
    "Color_Name",
    "Box_Color",
    "Median_Color",
    "Palette_Mode",
    "Hatch",
    "Box_Edge_Color",
)
MARKER_COLUMNS = (
    "Comparison_Block",
    "Group",
    "Well",
    "Technical_Well_Index",
    "Marker",
    "X_Offset",
    "Point_Color",
    "Point_Edge_Color",
)
HATCHES = ("", "//", "xx", "..", "\\\\", "++", "oo", "--", "**", "/o", "x.")


def prepare_palette(data: dict, *, survival: bool = False) -> dict:
    """Keep block-local roles and well identities stable across every view."""
    conditions, wells, palette_rows, marker_rows = {}, {}, [], []
    for block in data["blocks"]:
        groups = block["groups"]
        styles = condition_palette(groups, block["control"])
        for index, group in enumerate(groups):
            style = styles[group]
            # All template conditions count, even without measurements.
            hatch = ""
            if survival and len(groups) >= 6 and index:
                patterns = HATCHES[1:]
                hatch = patterns[(index - 1) % len(patterns)] * (
                    1 + (index - 1) // len(patterns)
                )
            style["hatch"] = hatch
            # Hatch strokes share the box edge color in Matplotlib. Use
            # a contrasting stroke on dark fills as well as light ones.
            style["edge_color"] = style["median_color"] if hatch else "#334155"
            palette_rows.append(
                {
                    "Comparison_Block": block["id"],
                    "Group": group,
                    "Is_Control": style["is_control"],
                    "Color_Name": style["color_name"],
                    "Box_Color": style["color"],
                    "Median_Color": style["median_color"],
                    "Palette_Mode": style["palette_mode"],
                    "Hatch": hatch,
                    "Box_Edge_Color": style["edge_color"],
                }
            )
            identities = sorted(
                row["Well"]
                for row in data["design"]
                if row["Comparison_Block"] == block["id"]
                and row["Group"] == group
            )
            for position, (well, marker) in enumerate(
                zip(identities, well_markers(len(identities)))
            ):
                record = {
                    "Comparison_Block": block["id"],
                    "Group": group,
                    "Well": well,
                    "Technical_Well_Index": position + 1,
                    "Marker": marker,
                    "X_Offset": (
                        -0.14 + 0.28 * position / (len(identities) - 1)
                        if len(identities) > 1
                        else 0.0
                    ),
                    "Point_Color": POINT_COLOR,
                    "Point_Edge_Color": POINT_EDGE_COLOR,
                }
                wells[well] = record
                marker_rows.append(record)
        conditions[block["id"]] = styles
    palette = {
        "conditions": conditions,
        "wells": wells,
        "palette_rows": palette_rows,
        "marker_rows": marker_rows,
        "survival_hatching_min_conditions": 6 if survival else None,
        "markers": well_markers(
            max((r["Technical_Well_Index"] for r in marker_rows), default=0)
        ),
    }
    data["plot_palette"] = palette
    return palette


def palette_tables(data):
    """Expose exact colors and full-template well assignments to Excel."""
    palette = data.get("plot_palette")
    if palette is None:
        return []
    return [
        ("Plot Palette", PALETTE_COLUMNS, palette["palette_rows"]),
        ("Well Markers", MARKER_COLUMNS, palette["marker_rows"]),
    ]


def save_palette(data, output, *, survival=False):
    palette = prepare_palette(data, survival=survival)
    save_json(output / "plot_palette.json", palette, allow_nan=False)


def draw_well_points(axis, records, metric, position, palette, *, size=48):
    """Draw measured wells only, without compacting gaps left by exclusions."""
    for row in records:
        style = palette["wells"][row["Well"]]
        axis.scatter(
            [position + style["X_Offset"]],
            [row[metric]],
            marker=style["Marker"],
            s=size,
            color=POINT_COLOR,
            edgecolors=POINT_EDGE_COLOR,
            linewidths=0.9,
            zorder=4,
            clip_on=False,
        )


def marker_legend_height(palette):
    return 0.4 + 0.3 * math.ceil(len(palette["markers"]) / 6)


def draw_marker_legend(axis, palette):
    """Explain ordinal shapes; the saved mapping identifies exact wells."""
    from matplotlib.lines import Line2D

    axis.set_axis_off()
    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker=marker,
            markersize=6,
            markerfacecolor=POINT_COLOR,
            markeredgecolor=POINT_EDGE_COLOR,
            label=f"Well {index}",
        )
        for index, marker in enumerate(palette["markers"], 1)
    ]
    if handles:
        axis.legend(
            handles=handles,
            loc="center",
            ncol=min(6, len(handles)),
            frameon=False,
            fontsize=8,
            title="Point shape: well within condition (see Well Markers)",
            title_fontsize=9,
        )
