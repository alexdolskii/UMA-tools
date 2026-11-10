"""Read immutable UMA outputs and resolve run-qualified display assets."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

from .files import sha256_file

PREFIX = "UMA_Web_Report_"
KINDS = {
    "Alignment_assay_results_": "Alignment",
    "Thickness_assay_results_": "Thickness",
    "Area_assay_results_": "Area",
    "Combined_Results_": "Collection",
    "UMA_Report_": "Report",
    "Report_Batch_": "Report plan",
}
COPY_SUFFIXES = {".csv", ".xlsx", ".json", ".pdf", ".png", ".log", ".txt"}


def uid(value):
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:24]


def basename(value):
    return str(value or "").replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]


def numeric(value):
    try:
        result = float(value)
    except (ValueError, TypeError):
        return None
    return result if math.isfinite(result) else None


def truth(value):
    return str(value).lower() in ("true", "1")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        if not columns or any(not x for x in columns):
            raise ValueError(f"Missing CSV headers: {path}")
        if len(set(columns)) != len(columns):
            raise ValueError(f"Duplicate CSV headers: {path}")
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f"CSV row/header length mismatch: {path}")
    return {"columns": columns, "rows": rows}


def regular(path, root):
    """Do not traverse links or AppleDouble/hidden entries."""
    try:
        relative = path.relative_to(root)
    except ValueError:
        return False
    current = root
    for part in relative.parts:
        current = current / part
        if part.startswith(".") or current.is_symlink():
            return False
    return path.is_file()


def walk_files(root):
    for child in sorted(root.iterdir()):
        if child.name.startswith(".") or child.is_symlink():
            continue
        if child.is_dir():
            yield from walk_files(child)
        elif child.is_file():
            yield child


def discover(source):
    """Only recognize analysis roots, never archived or exported copies."""
    assay = source / "uma_assay"
    if assay.is_symlink() or not assay.is_dir():
        raise ValueError(f"No regular uma_assay directory: {source}")
    result = []
    for path in sorted(assay.iterdir()):
        if not path.is_dir() or path.is_symlink():
            continue
        for prefix, kind in KINDS.items():
            if path.name.startswith(prefix):
                result.append((path, kind))
                break
    return result


class Resolver:
    """Relocate exact experiment-relative paths without basename searches."""

    def __init__(self, source):
        self.source = source.resolve()

    def resolve(self, value, base=None, directory=False):
        text = str(value or "").replace("\\", "/")
        if not text or ".." in text.split("/"):
            return None
        path = Path(text)
        candidates = []
        if path.is_absolute():
            candidates.append(path)
        elif not re.match(r"^[a-zA-Z]:/", text):
            candidates.append((base or self.source) / path)
        parts = text.split("/")
        for anchor in ("uma_assay", "fia_assay"):
            if parts.count(anchor) == 1:
                candidates.append(
                    self.source.joinpath(*parts[parts.index(anchor) :])
                )
        found = set()
        for candidate in candidates:
            if not candidate.is_relative_to(self.source):
                continue
            check = candidate / "dummy" if directory else candidate
            if any(
                p.startswith(".") for p in check.relative_to(self.source).parts
            ):
                continue
            current = self.source
            linked = False
            for part in candidate.relative_to(self.source).parts:
                current = current / part
                linked |= current.is_symlink()
            if linked:
                continue
            if candidate.is_dir() if directory else candidate.is_file():
                found.add(candidate.resolve())
        return next(iter(found)) if len(found) == 1 else None


def rows(run, name):
    return run["tables"].get(name.lower(), {}).get("rows", [])


def load_run(path, kind, output, logger):
    """Copy source documents byte-for-byte and expose exact CSV strings."""
    identity = uid(path.name)
    run = {
        "id": identity,
        "name": path.name,
        "kind": kind,
        "source": str(path),
        "tables": {},
        "documents": {},
        "downloads": [],
        "warnings": [],
        "hashes": {},
    }
    for file in walk_files(path):
        if file.suffix.lower() not in COPY_SUFFIXES:
            continue
        relative = file.relative_to(path).as_posix()
        try:
            before = sha256_file(file)
            content = file.read_bytes()
            if hashlib.sha256(content).hexdigest() != before:
                raise ValueError("Source changed while being copied")
            target = output / "downloads" / identity / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
            run["hashes"][relative] = before
            run["downloads"].append(
                {
                    "name": relative,
                    "path": target.relative_to(output).as_posix(),
                    "source": str(file),
                    "sha256": before,
                    "bytes": len(content),
                }
            )
            if file.suffix.lower() == ".csv" and relative.split("/")[
                0
            ] not in {"Inputs", "Raw"}:
                run["tables"][relative[:-4].lower()] = read_csv(target)
            elif file.suffix.lower() == ".json":
                run["documents"][relative] = read_json(target)
            if sha256_file(file) != before:
                raise ValueError("Source changed during export")
        except (OSError, ValueError, csv.Error) as error:
            run["warnings"].append(f"{relative}: {error}")
            logger.exception("DOCUMENT_FAILED | %s", file)
    status = run["documents"].get("run_status.json", {})
    run["status"] = status.get("status", "UNVERIFIED")
    run["parameters"] = run["documents"].get("run_parameters.json", {})
    run["route"] = (
        run["parameters"].get("fiber_route")
        or status.get("fiber_route")
        or run["documents"].get("plate_design.json", {}).get("Fiber_route")
        or {}
    )
    artifact_hash = status.get("artifacts_sha256")
    if artifact_hash:
        if run["hashes"].get("artifacts.json") != artifact_hash:
            run["warnings"].append("Recorded artifacts.json checksum mismatch")
        for record in run["documents"].get("artifacts.json", []):
            relative = record["path"].replace("\\", "/")
            if run["hashes"].get(relative) != record["sha256"]:
                run["warnings"].append(
                    f"Recorded checksum mismatch: {record['path']}"
                )
    for record in rows(run, "input_manifest"):
        relative = record.get("Archive_Relative_Path", "").replace("\\", "/")
        if relative and run["hashes"].get(relative) != record.get("SHA256"):
            run["warnings"].append(f"Archived checksum mismatch: {relative}")
    run["integrity_ok"] = not run["warnings"]
    return run


def field(row, name):
    return row.get("FN_Source__" + name, row.get(name, ""))


def observations(run):
    """Retain excluded and unmatched rows, without inventing zero values."""
    if run["kind"] == "Report":
        values = rows(run, "fia_image_audit") or rows(run, "merged_data")
    elif run["kind"] == "Area":
        values = (
            rows(run, "fibronectin_area_summary")
            or rows(run, "fibronectin_area_summary.partial")
            or rows(run, "fn_projection_manifest")
            or rows(run, "fn_projection_manifest.partial")
        )
    else:
        values = rows(run, run["kind"].lower() + "_summary")
        if run["kind"] == "Alignment":
            values = values or rows(run, "analysis/alignment_summary")
    result = []
    for row in values:
        name = row.get("Image_ID") or row.get("Source_File_Name")
        if not name:
            name = row.get("File_Name", "")
            if run["kind"] == "Alignment":
                # Historical alignment names omit the original extension;
                # show them without guessing an ND2/TIFF source identity.
                name = name.removesuffix(
                    "_processed_orientation_distribution.csv"
                )
        if not name:
            continue
        result.append(
            {
                "image": name,
                "row": row,
                "group": row.get("Group", ""),
                "well": row.get("Well", ""),
                "files": [],
                "warnings": [],
                "included": truth(row.get("Included_In_Filtered_Plots")),
                "reason": row.get("Exclusion_Reason", ""),
                "source": field(row, "Source_Original_Path"),
            }
        )
    return result


def link_observation(run, observation, assets):
    """Use the collection's frozen selection, never the latest measurement."""
    base = Path(run["source"])
    row, name = observation["row"], observation["image"]
    original_stem = Path(name).stem if Path(name).suffix else name
    route = run["route"]
    cohort = json.dumps(route, sort_keys=True)
    selected = (
        run["documents"]
        .get("Inputs/collector_run_status.json", {})
        .get("selected", {})
    )

    def add(value, parent, label, mode="linear", expected=None, group=""):
        observation["files"].append(
            assets.add(
                value,
                parent,
                label,
                f"{parent}|{cohort}|{group or label}",
                mode,
                expected,
            )
        )

    if run["kind"] in ("Report", "Area"):
        area_base = (
            assets.resolver.resolve(
                selected.get("Area", {}).get("run_directory"), directory=True
            )
            if run["kind"] == "Report"
            else base
        )
        for column, label, mode, hash_name in (
            (
                "Source_Projection_Path",
                "Area SUM projection",
                "linear",
                "Projection_SHA256",
            ),
            ("Mask_Path", "Area mask", "binary", "Mask_SHA256"),
        ):
            value = field(row, column)
            if value:
                add(
                    value,
                    area_base or base,
                    label,
                    mode,
                    field(row, hash_name) or None,
                    label + "|Z=" + str(field(row, "Number_of_Z_Stacks")),
                )
        digest = field(row, "Source_SHA256")
        if digest and route.get("marker_id") and route.get("source_channel"):
            observation["comparison_key"] = uid(
                f"{name}|{digest}|{route['marker_id']}|"
                f"{route['source_channel']}"
            )
        observation["threshold"] = {
            key: field(row, key)
            for key in (
                "Threshold_Lower",
                "Threshold_Upper",
                "Effective_Threshold_Upper",
                "Threshold_Units",
            )
        }
    for kind in ("Alignment", "Thickness"):
        parent = (
            base
            if run["kind"] == kind
            else assets.resolver.resolve(
                selected.get(kind, {}).get("run_directory"), directory=True
            )
        )
        if parent is None:
            if run["kind"] == "Report":
                observation["warnings"].append(
                    f"{kind}: recorded measurement run unavailable"
                )
            continue
        if run["kind"] not in (kind, "Report"):
            continue
        if kind == "Alignment":
            stem = original_stem + "_processed"
            add(
                f"{stem}.tif",
                parent,
                "Alignment MAX projection",
                group="Alignment MAX projection|Z="
                + str(row.get("Number_of_Z_Stacks", "unrecorded")),
            )
            add(
                f"Images/{stem}_orientation_composition.png",
                parent,
                "Fiber orientation",
                "color",
            )
            add(
                f"Images/normalized_images/{stem}_normalized_orientation.png",
                parent,
                "Normalized fiber orientation",
                "color",
            )
        else:
            add(f"Mask_{name}.tif", parent, "Thickness XZ mask", "binary")
            add(f"Local_Thickness_{name}.tif", parent, "Local thickness map")
            observation["warnings"].append(
                "Unmasked XZ projection is not saved by this workflow."
            )
    link_nuclei(run, observation, assets, add)


def link_nuclei(run, observation, assets, add):
    """Show only the nuclei run already matched by UMA, with proven lineage."""
    row = observation["row"]
    if row.get("FIA_Match_Status") != "MATCHED":
        return
    selection = run["documents"].get("fia_selection.json", {})
    parent = assets.resolver.resolve(selection.get("run"), directory=True)
    mask = row.get("FIA_Mask_Name", "")
    if parent is None or parent.name != row.get("FIA_Run_ID"):
        observation["warnings"].append("Matched FIA run unavailable")
        return
    snapshots = [
        content
        for name, content in run["documents"].items()
        if name.startswith("Inputs/FIA/") and name.endswith("_nuclei_run.json")
    ]
    current = assets.resolver.resolve("nuclei_run.json", parent)
    if len(snapshots) != 1 or current is None:
        observation["warnings"].append("FIA run provenance unavailable")
        return
    try:
        if read_json(current) != snapshots[0]:
            raise ValueError("FIA metadata changed since the UMA report")
        if not mask or basename(mask) != mask:
            raise ValueError("Invalid recorded FIA mask name")
        add(mask, parent, "FIA nuclear mask", "binary")
        star = assets.resolver.resolve(
            snapshots[0].get("stardist_folder"), directory=True
        )
        record = (
            assets.resolver.resolve("stardist_run.json", star)
            if (star is not None)
            else None
        )
        image = mask.removesuffix("_StarDist_processed_processed.tif") + ".tif"
        checksum = (
            read_json(record).get("sources", {}).get(image)
            if (record is not None)
            else None
        )
        if not isinstance(checksum, str) or len(checksum) != 64:
            raise ValueError("Nuclei projection has no verifiable source hash")
        projection = assets.resolver.source / "fia_assay" / "Nuclei" / image
        add(
            str(projection), parent, "FIA nuclei projection", expected=checksum
        )
    except (OSError, ValueError, TypeError) as error:
        observation["warnings"].append(str(error))
