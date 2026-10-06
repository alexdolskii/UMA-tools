"""Validate optional plate-map display order independently of cell roles."""

from __future__ import annotations

import math
import re

from .report_schema import ValidationError

GROUP_ORDER_COLUMNS = (
    "Group",
    "Order",
    "Source",
    "Order_Cell",
    "Group_Cell",
    "Grid_Cells",
)


def _invalid(sheet, location, issue, group=""):
    address = f"{sheet.title}!{location}"
    raise ValidationError(
        f"{address}: {issue}",
        [
            {
                "Table": "Group_Order",
                "Excel_Cell": address,
                "Group": group,
                "Issue": issue,
            }
        ],
    )


def _positive_integer(cell):
    value = cell.value
    if cell.data_type in ("f", "e") or isinstance(value, bool):
        return None
    if isinstance(value, str):
        if not re.fullmatch(r"[0-9]+", value.strip()):
            return None
        value = int(value.strip())
    if isinstance(value, int):
        return value if value > 0 else None
    if isinstance(value, float) and math.isfinite(value):
        return int(value) if value > 0 and value.is_integer() else None
    return None


def read_group_order(sheet, well_map, cells):
    """Read ranks to the right of A1:M9, or retain row-major grid order.

    A condition name appears once in the order table even if it belongs
    to several functional comparison blocks. Formatting in this table
    never establishes block membership or a control role.
    """
    groups = list(dict.fromkeys(well_map.values()))
    columns = {"order": [], "group": []}
    for row in sheet.iter_rows(min_row=1, max_row=1, min_col=14):
        for cell in row:
            header = str(cell.value or "").strip().casefold()
            if header in {"order", "group", "groups"}:
                columns["order" if header == "order" else "group"].append(cell)
    explicit, numbers = {}, {}
    if any(columns.values()):
        if any(len(found) != 1 for found in columns.values()):
            addresses = ", ".join(
                c.coordinate for found in columns.values() for c in found
            )
            _invalid(
                sheet,
                addresses,
                "expected one Order and one Group (or Groups) header "
                "in row 1 to the right of the plate grid",
            )
        order_col, group_col = (
            columns[key][0].column for key in ("order", "group")
        )
        for merged in sheet.merged_cells.ranges:
            if any(
                merged.min_col <= c <= merged.max_col
                for c in (order_col, group_col)
            ):
                _invalid(
                    sheet,
                    str(merged),
                    "do not merge cells in Order/Group columns",
                )
        for row in range(2, sheet.max_row + 1):
            number, name = (
                sheet.cell(row, order_col),
                sheet.cell(row, group_col),
            )
            blank = [
                c.value is None
                or (isinstance(c.value, str) and not c.value.strip())
                for c in (number, name)
            ]
            if all(blank):
                continue
            if any(blank):
                _invalid(
                    sheet,
                    f"{number.coordinate}/{name.coordinate}",
                    "fill both Order and Group, or leave both blank",
                )
            if name.data_type in ("f", "e") or not isinstance(name.value, str):
                _invalid(
                    sheet,
                    name.coordinate,
                    "Group must be literal text matching the plate grid",
                )
            rank = _positive_integer(number)
            if rank is None:
                _invalid(
                    sheet,
                    number.coordinate,
                    "Order must be a positive whole number",
                    name.value,
                )
            if name.value not in groups:
                _invalid(
                    sheet,
                    name.coordinate,
                    f"unknown Group {name.value!r}; use the exact "
                    "condition name from the plate grid",
                    name.value,
                )
            if name.value in explicit:
                previous = explicit[name.value]["Group_Cell"]
                _invalid(
                    sheet,
                    name.coordinate,
                    f"repeated Group {name.value!r}; already listed at "
                    f"{previous}",
                    name.value,
                )
            if rank in numbers:
                _invalid(
                    sheet,
                    number.coordinate,
                    f"repeated Order {rank}; already used at {numbers[rank]}",
                    name.value,
                )
            numbers[rank] = f"{sheet.title}!{number.coordinate}"
            explicit[name.value] = {
                "Group": name.value,
                "Order": rank,
                "Source": "order_table",
                "Order_Cell": f"{sheet.title}!{number.coordinate}",
                "Group_Cell": f"{sheet.title}!{name.coordinate}",
            }
        missing = [g for g in groups if g not in explicit]
        if explicit and missing:
            _invalid(
                sheet,
                "Order/Group",
                "incomplete order table; add conditions: "
                + ", ".join(map(repr, missing)),
            )
    records = []
    for index, group in enumerate(groups, 1):
        record = explicit.get(
            group,
            {
                "Group": group,
                "Order": index,
                "Source": "plate_grid",
                "Order_Cell": "",
                "Group_Cell": "",
            },
        )
        records.append(
            {
                **record,
                "Grid_Cells": "; ".join(
                    f"{sheet.title}!{cells[well]}"
                    for well, label in well_map.items()
                    if label == group
                ),
            }
        )
    return sorted(records, key=lambda record: record["Order"])


def order_blocks(blocks, records):
    """Reorder presentation while preserving block IDs and control roles."""
    ranks = {row["Group"]: row["Order"] for row in records}
    ordered = [
        {**block, "groups": sorted(block["groups"], key=ranks.__getitem__)}
        for block in blocks
    ]
    return sorted(
        ordered, key=lambda block: min(ranks[g] for g in block["groups"])
    )
