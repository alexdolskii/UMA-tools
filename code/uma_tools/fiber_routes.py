"""Keep one marker and original channel together throughout the UMA workflow.

The context is scoped to one sequential route, like ImageRun's image context.
It never changes scientific parameters or the physical source channel numbers.
"""

from __future__ import annotations

import json
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path

from . import plate_templates
from .contracts import IMAGE_EXTENSIONS
from .files import save_csv, save_json, sha256_file
from .plate_design import PlateDesign, digest, read_design, well_from_name
from .progress import console

_ACTIVE = ContextVar("uma_fiber_route", default=None)
SCHEMA = "uma-fiber-route-v1"


@dataclass(frozen=True)
class FiberRoute:
    source: Path
    design: PlateDesign
    marker: dict
    channel: int

    @property
    def wells(self):
        return {
            well: item
            for well, item in self.design.wells.items()
            if self.design.assigned(well, "Fiber", self.marker["id"])
            and self.design.channel(well, "Fiber", self.marker["id"])
            == self.channel
        }

    @property
    def label(self):
        # Marker IDs include a hash, so punctuation/Unicode cannot collide.
        return f"{self.marker['id']}_C{self.channel}"

    @property
    def metadata(self):
        return {
            "schema": SCHEMA,
            "marker_id": self.marker["id"],
            "marker_name": self.marker["name"],
            "method": "Fiber",
            "source_channel": self.channel,
            "scope_sha256": digest(
                {
                    well: {
                        "channel": self.channel,
                        "channel_count": sum(
                            c is not None for c in item["channels"]
                        ),
                    }
                    for well, item in self.wells.items()
                }
            ),
        }

    def verify(self):
        if sha256_file(self.design.path) != self.design.workbook_sha256:
            raise ValueError(
                f"Plate scheme changed during command: {self.design.path}"
            )

    def inventory(self):
        records, stems = [], set()
        for path in sorted(self.source.iterdir()):
            if (
                path.name.startswith(".")
                or path.is_symlink()
                or not path.is_file()
                or path.suffix.lower() not in IMAGE_EXTENSIONS
            ):
                continue
            well = well_from_name(path.name)
            selected = well in self.wells
            if selected and path.stem in stems:
                raise ValueError(f"Ambiguous original image stem: {path.stem}")
            if selected:
                stems.add(path.stem)
            records.append(
                {
                    "File_Name": path.name,
                    "Well": well,
                    "Marker": self.marker["name"],
                    "Marker_ID": self.marker["id"],
                    "Source_Channel": self.channel if selected else "",
                    "Status": "SELECTED" if selected else "NOT_APPLICABLE",
                    "Reason": "Assigned marker and channel"
                    if selected
                    else (
                        "Well is unannotated or this marker/channel "
                        "is not assigned"
                    ),
                }
            )
        return records


def routes(source, design, requested=None):
    markers = design.markers("Fiber")
    if requested:
        markers = [m for m in markers if requested in (m["id"], m["name"])]
    return [
        FiberRoute(Path(source).resolve(), design, marker, channel)
        for marker in markers
        for channel in sorted(set(marker["channels"].values()))
    ]


def preflight(folders, args, validator=None):
    """Validate every folder and source identity before starting any engine."""
    designs = {}

    def validate(path, source):
        design = read_design(path)
        if not design.markers("Fiber"):
            raise ValueError("No Fiber: marker name assignments in Channels")
        # Check names and duplicate stems without loading microscopy pixels.
        for route in routes(source, design):
            route.inventory()
        if validator:
            validator(path, source)
        designs[source.resolve()] = design
        return "Plate Map"

    prepared = plate_templates.preflight_templates(
        [Path(p).resolve() for p in folders],
        validate,
        interactive=plate_templates.interactive_selection(args),
        prefix="FIA",
    )
    requested = getattr(args, "marker", None)
    if requested and not any(
        routes(source, design, requested) for source, design in designs.items()
    ):
        raise ValueError(
            f"Fiber marker {requested!r} is not declared in these schemes"
        )
    return [(item, designs.get(item.source.resolve())) for item in prepared]


def active_route():
    return _ACTIVE.get()


@contextmanager
def use_route(route):
    route.verify()
    token = _ACTIVE.set(route)
    console(
        f"Fiber: {route.marker['name']} | original channel C{route.channel} | "
        f"{route.source.name}"
    )
    try:
        yield route
        route.verify()
    finally:
        _ACTIVE.reset(token)


def selected_name(source, name):
    route = active_route()
    return (
        route is None
        or Path(source).resolve() != route.source
        or well_from_name(name) in route.wells
    )


def validate_channels(name, count):
    route = active_route()
    if route is not None:
        well = well_from_name(name)
        if well not in route.wells:
            raise ValueError(
                f"Image does not belong to selected Fiber route: {name}"
            )
        route.design.validate_channels(well, count)


def snapshot(output, state=None):
    route = active_route()
    if route is None:
        return
    route.verify()
    output = Path(output)
    (output / "plate_design.xlsx").write_bytes(route.design.content)
    save_json(
        output / "plate_design.json",
        {
            "Schema_version": 1,
            "Workbook": str(route.design.path),
            "Workbook_SHA256": route.design.workbook_sha256,
            "Wells": route.design.wells,
            "Panels": route.design.panels,
            "Nuclei_SHA256": route.design.scope_hash("Nuclei"),
            "Fiber_route": route.metadata,
        },
    )
    records = route.inventory()
    if records:
        save_csv(output / "source_selection.csv", list(records[0]), records)
    if state is not None:
        state["fiber_route"] = route.metadata


def validate_run(folder, state):
    """Reject old/mismatched measurements; formatting edits remain reusable."""
    route = active_route()
    if route is None:
        return
    if state.get("fiber_route") != route.metadata:
        raise ValueError(
            "Missing or incompatible Fiber marker/channel provenance; rerun "
            "the affected analysis"
        )
    path = Path(folder)
    meta = json.loads((path / "plate_design.json").read_text())
    original = read_design(path / "plate_design.xlsx")
    if (
        meta.get("Schema_version") != 1
        or meta.get("Workbook_SHA256") != original.workbook_sha256
        or meta.get("Wells") != original.wells
        or meta.get("Panels") != original.panels
        or meta.get("Fiber_route") != route.metadata
    ):
        raise ValueError(f"Invalid plate scheme snapshot: {folder}")
    previous = routes(route.source, original, route.marker["id"])
    if not any(r.metadata == route.metadata for r in previous):
        raise ValueError(f"Fiber assignments have changed: {folder}")


def annotate_summary(path, source, records):
    """Attach source identity without altering measured fields."""
    import csv

    route = active_route()
    if route is None:
        return
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns, rows = list(reader.fieldnames), list(reader)
    suffix = "_processed_orientation_distribution.csv"
    lookup = {
        Path(r["File_Name"]).stem: r["File_Name"]
        for r in route.inventory()
        if r["Status"] == "SELECTED"
    }
    for row in rows:
        name = row["File_Name"]
        if name.endswith(suffix):
            name = lookup[name.removesuffix(suffix)]
        original = Path(source) / name
        source_hash = sha256_file(original)
        if records[name].get("Source_SHA256") != source_hash:
            raise ValueError(
                f"Original changed during measurement: {original}"
            )
        if row.get("Source_SHA256") and row["Source_SHA256"] != source_hash:
            raise ValueError(
                f"Original changed during measurement: {original}"
            )
        row.update(
            Fiber_Marker=route.marker["name"],
            Fiber_Marker_ID=route.marker["id"],
            Source_Channel=route.channel,
            Source_File_Name=name,
            Source_SHA256=source_hash,
        )
    extra = [
        "Fiber_Marker",
        "Fiber_Marker_ID",
        "Source_Channel",
        "Source_File_Name",
        "Source_SHA256",
    ]
    save_csv(
        Path(path), columns + [c for c in extra if c not in columns], rows
    )


def report_metadata(data):
    route = active_route()
    if route is not None:
        data.update(
            fiber_route=route.metadata,
            layout_panels=route.design.panels,
            applicable_wells=list(route.wells),
        )
