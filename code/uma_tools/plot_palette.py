"""Box colors, neutral image points and well markers for reports.

This is a new selection of historical digital colors, not a numbered
Wada combination. Green tints are derived adaptations. See the README
for the palette order and sources.
"""

CONTROL_COLOR = ("Grayish Lavender A", "#B5B1D8")
TREATMENT_COLORS = (
    ("Dusky Green", "#004F46"),
    ("Orange", "#F37420"),
    ("Deep Indigo", "#051230"),
    ("Dull Blue Violet", "#80719E"),
    ("Ivory Buff", "#EBD3A2"),
    ("Violet", "#4F4086"),
    ("Verditter Blue", "#6FB5A8"),
)
POINT_COLOR = "#D0D0D0"
POINT_EDGE_COLOR = "#111314"
LOW_FN_EDGE_COLOR = "#D62728"
WELL_MARKERS = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "p", "h", "8")


def median_color(fill):
    """Choose a light or dark median to contrast with the box fill."""

    def luminance(color):
        rgb = [int(color[i : i + 2], 16) / 255 for i in (1, 3, 5)]
        linear = [
            value / 12.92
            if value <= 0.04045
            else ((value + 0.055) / 1.055) ** 2.4
            for value in rgb
        ]
        return sum(v * w for v, w in zip(linear, (0.2126, 0.7152, 0.0722)))

    background, dark = luminance(fill), luminance(POINT_EDGE_COLOR)
    light_contrast = 1.05 / (background + 0.05)
    dark_contrast = (max(background, dark) + 0.05) / (
        min(background, dark) + 0.05
    )
    return "#FFFFFF" if light_contrast > dark_contrast else POINT_EDGE_COLOR


def well_markers(count):
    """Use shapes, then numbers; never cycle within a condition."""
    return [
        WELL_MARKERS[index] if index < len(WELL_MARKERS) else f"${index + 1}$"
        for index in range(count)
    ]


def _green_tints(count):
    """Blend Dusky Green toward white in sRGB, capped at 75% white."""
    base_name, base_hex = TREATMENT_COLORS[0]
    base = tuple(int(base_hex[i : i + 2], 16) for i in (1, 3, 5))
    colors = []
    for index in range(count):
        fraction = 0.75 * index / max(1, count - 1)
        color = "#" + "".join(
            f"{round(value + (255 - value) * fraction):02X}" for value in base
        )
        name = base_name if index == 0 else f"{base_name} tint {index + 1}"
        colors.append((name, color))
    return colors


def condition_palette(groups, control=None):
    """Reserve lavender for a known control; keep condition order."""
    treatments = [group for group in groups if group != control]
    # Without a valid control, lavender remains unused. Eight unmarked
    # conditions therefore need tints rather than a false control color.
    tinted = len(groups) > 8 or len(treatments) > len(TREATMENT_COLORS)
    colors = _green_tints(len(treatments)) if tinted else TREATMENT_COLORS
    assigned = dict(zip(treatments, colors))
    if control is not None:
        assigned[control] = CONTROL_COLOR
    return {
        group: {
            "color_name": assigned[group][0],
            "color": assigned[group][1],
            "median_color": median_color(assigned[group][1]),
            "is_control": group == control,
            "palette_mode": "green_tints" if tinted else "historical",
        }
        for group in groups
    }


def _block_control(groups, design, expected_wells):
    """Accept one bold condition with consistent whole-cell styles."""
    controls = []
    for group in groups:
        records = [row for row in design if row["Group"] == group]
        identities = {
            (row.get("Color_Code"), row.get("Is_Control")) for row in records
        }
        if {row["Well"] for row in records} != expected_wells[group] or len(
            identities
        ) != 1:
            return None
        _, bold = identities.pop()
        if bold is None:
            return None
        if bold:
            controls.append(group)
    return controls[0] if len(controls) == 1 else None


def report_palette(data, panels, log):
    """Assign styles before filtering, including unimaged controls."""
    statistics = data.get("statistics") or {}
    design = data.get("plot_design", statistics.get("design", []))
    well_map = data.get("well_map") or {
        well: group
        for group, wells in data["group_wells"].items()
        for well in wells
    }
    expected_wells = {
        group: {well for well, label in well_map.items() if label == group}
        for group in well_map.values()
    }
    conditions, blocks = {}, []
    for panel in panels:
        groups = list(panel["groups"])
        # Count all template conditions, even when one has
        # no source images. The control reserves its lavender slot.
        for row in design:
            if (
                row.get("Color_Code") == panel["color_code"]
                and row["Group"] not in groups
            ):
                groups.append(row["Group"])
        control = _block_control(groups, design, expected_wells)
        styles = condition_palette(groups, control)
        conditions.update(styles)
        blocks.append(
            {
                "panel": panel["id"],
                "color_code": panel["color_code"],
                "groups": groups,
                "control": control,
                "palette_mode": styles[groups[0]]["palette_mode"],
            }
        )
        if control is None:
            log.event(
                "WARNING",
                "Plot palette",
                f"Panel {panel['id']}: no unambiguous bold control. "
                "Lavender is reserved; colors describe conditions only.",
            )
    markers = well_markers(max(map(len, data["group_wells"].values())))
    return {
        "point_color": POINT_COLOR,
        "point_edge_color": POINT_EDGE_COLOR,
        "conditions": conditions,
        "blocks": blocks,
        "markers": markers,
        "wells": {
            well: markers[index]
            for wells in data["group_wells"].values()
            for index, well in enumerate(wells)
        },
    }
