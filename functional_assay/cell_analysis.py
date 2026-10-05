"""Analyze stitched wells and retain masks, measurements, and provenance."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections.abc import Callable, Sequence
from importlib.metadata import version
from pathlib import Path
from typing import Any

from uma_tools.cli import _run_imagej_command
from uma_tools.config import read_config
from uma_tools.files import safe_label, save_csv, save_json, sha256_file
from uma_tools.imagej import FIJI_ENDPOINT, initialize_imagej
from uma_tools.run import unique_output, utc_now

from .calibration import read_well_calibration
from .stitching import discover_wells, validate_frames
from .workflow import (
    NoInputError,
    RunLog,
    activity,
    assay_directory,
    batch_error,
    command_error,
    exclusions,
    outcome,
    save_exclusions,
)

SUMMARY_COLUMNS = (
    "Well",
    "File_Name",
    "Status",
    "Error",
    "Object_Count",
    "Mask_Area_px2",
    "Mask_Area_um2",
    "Counted_Object_Area_px2",
    "Counted_Object_Area_um2",
    "Threshold_Method",
    "Threshold_Lower",
    "Threshold_Upper",
    "Min_Size_px2",
    "Min_Size_um2",
    "Pixel_Size_X_um",
    "Pixel_Size_Y_um",
    "Width_px",
    "Height_px",
    "Width_um",
    "Height_um",
    "Overlap_Percent",
    "Overlap_Status",
)
OBJECT_COLUMNS = ("Well", "Object_ID", "Area_px2", "Area_um2")
STITCHED_PATTERN = re.compile(r"(Well[A-Za-z0-9]+)_stitched\.tif", re.I)


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Validate command options before loading Fiji or scientific modules."""
    parser = argparse.ArgumentParser(
        description="Measure mask area and count cells in stitched wells"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('uma-functional-assay')}",
    )
    parser.add_argument("-i", "--input", required=True, help="UMA input JSON")
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        nargs="*",
        metavar="VALUE",
        help=(
            "Omit for RenyiEntropy. Bare -t uses 50 65535; "
            "-t LOWER [UPPER] sets manual 16-bit thresholds"
        ),
    )
    size = parser.add_mutually_exclusive_group()
    size.add_argument(
        "--min-size-px",
        type=float,
        help="Minimum particle area in pixels squared (default: 5)",
    )
    size.add_argument(
        "--min-size-um2",
        type=float,
        help="Minimum particle area in square micrometers",
    )
    args = parser.parse_args(argv)
    if args.threshold is not None:
        if len(args.threshold) > 2:
            parser.error("-t accepts at most LOWER and UPPER")
        bounds = (args.threshold or [50]) + [65535]
        lower, upper = bounds[:2]
        if not (
            math.isfinite(lower)
            and math.isfinite(upper)
            and 0 <= lower <= upper <= 65535
        ):
            parser.error(
                "Thresholds must satisfy 0 <= LOWER <= UPPER <= 65535"
            )
        args.threshold = (lower, upper)
    if args.min_size_px is None and args.min_size_um2 is None:
        args.min_size_px = 5.0
    for name in ("min_size_px", "min_size_um2"):
        value = getattr(args, name)
        if value is not None and (not math.isfinite(value) or value < 0):
            parser.error("Minimum particle area must be finite and >= 0")
    return args


def minimum_sizes(args: argparse.Namespace, pixel_area: float) -> tuple:
    """Return the same area cutoff in both pixel and physical units."""
    if not math.isfinite(pixel_area) or pixel_area <= 0:
        raise ValueError("Physical pixel area must be finite and positive")
    if args.min_size_um2 is not None:
        minimum_px = args.min_size_um2 / pixel_area
        # Avoid rejecting an exactly N-pixel object through unit rounding.
        if math.isclose(
            minimum_px, round(minimum_px), abs_tol=1e-12, rel_tol=0
        ):
            minimum_px = float(round(minimum_px))
        return minimum_px, args.min_size_um2
    return args.min_size_px, args.min_size_px * pixel_area


def discover_stitched(folder: Path) -> dict[str, Path]:
    """Find only visible stitched TIFFs; reject ambiguous duplicate wells."""
    if folder.is_symlink() or not folder.is_dir():
        raise NoInputError(f"Stitching results not found: {folder}")
    images = {}
    for path in sorted(folder.iterdir()):
        if path.name.startswith(".") or not path.is_file():
            continue
        match = STITCHED_PATTERN.fullmatch(path.name)
        if match is None:
            continue
        well = match[1]
        if any(key.casefold() == well.casefold() for key in images):
            raise ValueError(f"Duplicate stitched TIFF for {well}")
        images[well] = path
    if not images:
        raise NoInputError(f"No *_stitched.tif images found in {folder}")
    return images


def read_stitching_metadata(folder: Path) -> dict:
    """Require a finalized, auditable new-layout stitching run."""
    path = folder / "stitching_metadata.json"
    if folder.is_symlink() or path.is_symlink() or not path.is_file():
        raise NoInputError(f"Missing regular stitching metadata: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or not isinstance(
        data.get("wells"), dict
    ):
        raise ValueError("Unsupported or invalid stitching_metadata.json")
    records = data["wells"]
    completed = sum(r.get("status") == "completed" for r in records.values())
    if data.get("status") not in {"SUCCESS", "PARTIAL"}:
        raise NoInputError(
            f"Stitching run is not usable: {data.get('status')}"
        )
    if (
        not completed
        or data.get("completed_wells") != completed
        or data.get("failures") != len(records) - completed
        or any(
            r.get("status") not in {"completed", "failed"}
            for r in records.values()
        )
        or data["status"] != outcome(completed, len(records) - completed)
    ):
        raise ValueError("Inconsistent stitching completion record")
    overlap = data.get("overlap_percent")
    if (
        not isinstance(overlap, (int, float))
        or not math.isfinite(overlap)
        or not 0 <= overlap < 100
    ):
        raise ValueError("Invalid overlap in stitching_metadata.json")
    return data


def verify_stitching_record(
    path: Path, well: str, metadata: dict | None, calibration: dict
) -> dict[str, Any]:
    """Bind recorded overlap to this TIFF, and recheck the original scales."""
    digest = sha256_file(path)
    if metadata is None:
        return {"sha256": digest, "overlap_percent": None}
    record = metadata["wells"].get(well)
    if not record or record.get("status") != "completed":
        raise ValueError(f"No successful stitching metadata for {well}")
    if record.get("sha256") != digest:
        raise ValueError(f"{well}: TIFF does not match stitching metadata")
    recorded_scale = record.get("calibration")
    if recorded_scale is not None:
        for axis in ("x", "y"):
            key = f"pixel_size_{axis}_um"
            if not math.isclose(
                recorded_scale[key], calibration[key], rel_tol=1e-6
            ):
                raise ValueError(f"{well}: original calibration has changed")
    recorded_names = sorted(
        item["filename"] for item in record["source_frames"]
    )
    actual_names = sorted(item["filename"] for item in calibration["tiles"])
    if recorded_names != actual_names:
        raise ValueError(f"{well}: original filenames differ from stitching")
    return {"sha256": digest, "overlap_percent": metadata["overlap_percent"]}


def measure_well(
    well: str,
    path: Path,
    files: list,
    output: Path,
    args: argparse.Namespace,
    metadata: dict | None,
) -> tuple[dict, dict]:
    """Validate one well, analyze it, then write its scientific outputs."""
    from .cell_imagej import analyze_image, save_contours, save_mask

    validate_frames(files)
    activity(f"{well}: original calibration and TIFF provenance")
    calibration = read_well_calibration(files)
    provenance = verify_stitching_record(path, well, metadata, calibration)
    minimum_px, minimum_um2 = minimum_sizes(
        args, calibration["pixel_area_um2"]
    )
    activity(f"{well}: projection, filters, masks and object detection")
    result = analyze_image(path, args.threshold, minimum_px)
    width, height = result["width_px"], result["height_px"]
    if metadata is not None:
        record = metadata["wells"][well]
        if (record["width_px"], record["height_px"]) != (width, height):
            raise ValueError(f"{well}: TIFF dimensions differ from metadata")
    pixel_area = calibration["pixel_area_um2"]
    areas = result["object_areas_px2"]
    mask_area = int(result["area_mask"].sum())
    counted_area = int(areas.sum())
    activity(f"{well}: save masks, contours and object areas")
    for name in ("area_mask", "counting_mask"):
        save_mask(output / f"{well}_{name}.tif", result[name], calibration)
    save_contours(output / f"{well}_counted_contours.png", result)
    save_csv(
        output / f"{well}_objects.csv",
        OBJECT_COLUMNS,
        (
            {
                "Well": well,
                "Object_ID": index,
                "Area_px2": int(area),
                "Area_um2": int(area) * pixel_area,
            }
            for index, area in enumerate(areas, 1)
        ),
    )
    row = {
        "Well": well,
        "File_Name": path.name,
        "Status": "completed",
        "Error": "",
        "Object_Count": len(areas),
        "Mask_Area_px2": mask_area,
        "Mask_Area_um2": mask_area * pixel_area,
        "Counted_Object_Area_px2": counted_area,
        "Counted_Object_Area_um2": counted_area * pixel_area,
        "Threshold_Method": "RenyiEntropy"
        if args.threshold is None
        else "manual",
        "Threshold_Lower": result["threshold_lower"],
        "Threshold_Upper": result["threshold_upper"],
        "Min_Size_px2": minimum_px,
        "Min_Size_um2": minimum_um2,
        "Pixel_Size_X_um": calibration["pixel_size_x_um"],
        "Pixel_Size_Y_um": calibration["pixel_size_y_um"],
        "Width_px": width,
        "Height_px": height,
        "Width_um": width * calibration["pixel_size_x_um"],
        "Height_um": height * calibration["pixel_size_y_um"],
        "Overlap_Percent": provenance["overlap_percent"],
        "Overlap_Status": "not_recorded" if metadata is None else "recorded",
    }
    return row, {"calibration": calibration, "stitched_tiff": provenance}


def process_folder(
    folder: Path, args: argparse.Namespace, ensure_context: Callable[[], None]
) -> tuple[int, int]:
    """Write a fresh run; preserve diagnostics and continue past a bad well."""
    _, output = unique_output(
        assay_directory(folder), f"Cell_Analysis_{safe_label(folder.name)}_"
    )
    log = RunLog(output, folder, "cell_count")
    rows, records = [], {}
    failures = 0
    status = {
        "status": "RUNNING",
        "started_utc": utc_now(),
        "source": str(folder),
        "output": str(output),
        "functional_assay_version": version("uma-functional-assay"),
        "uma_tools_version": version("uma-tools"),
        "fiji_endpoint": FIJI_ENDPOINT,
        "parameters": vars(args),
        "wells": records,
    }
    status_path = output / "run_status.json"

    def checkpoint():
        summary = output / "Cell_Analysis_Summary.csv"
        save_csv(summary, SUMMARY_COLUMNS, rows)
        status["summary_sha256"] = sha256_file(summary)
        save_json(status_path, status, allow_nan=False)

    try:
        save_json(status_path, status, allow_nan=False)
        log.event("STARTED", "Cell analysis", folder)
        log.event("INFO", "Output", output)
        log.event(
            "INFO", "Plan", "1/3 inputs; 2/3 measure wells; 3/3 save audit"
        )
        log.phase(1, 3, "Validate stitching")
        stitched = assay_directory(folder) / "Stitched_Results"
        metadata = read_stitching_metadata(stitched)
        records.update(
            {well: {"status": "pending"} for well in metadata["wells"]}
        )
        shutil.copy2(stitched / "stitching_metadata.json", output)
        images = discover_stitched(stitched)
        unknown = set(images) - metadata["wells"].keys()
        if unknown:
            raise ValueError(
                f"TIFFs absent from stitching audit: {sorted(unknown)}"
            )
        status["stitching_status"] = metadata["status"]
        originals = discover_wells(folder)
        for index, (well, stitched_record) in enumerate(
            metadata["wells"].items()
        ):
            path = images.get(well, stitched / f"{well}_stitched.tif")
            log.phase(
                2,
                3,
                "Measure",
                finished=index,
                count=len(metadata["wells"]),
                detail=well,
            )
            try:
                if stitched_record["status"] != "completed":
                    raise ValueError(
                        "Excluded by stitching: "
                        + stitched_record.get("error", "Failed well")
                    )
                activity(f"{well}: initialize ImageJ if needed")
                ensure_context()
                files = [
                    item
                    for name, items in originals.items()
                    if name.casefold() == well.casefold()
                    for item in items
                ]
                row, record = measure_well(
                    well, path, files, output, args, metadata
                )
                rows.append(row)
                records[well] = {"status": "completed", **record}
                log.event(
                    "INFO",
                    well,
                    (
                        f"{row['Object_Count']} objects; "
                        f"mask {row['Mask_Area_um2']} µm²; "
                        f"threshold {row['Threshold_Lower']}–"
                        f"{row['Threshold_Upper']}"
                    ),
                    console=False,
                )
            except Exception as error:
                failures += 1
                message = str(error) or type(error).__name__
                rows.append(
                    {
                        "Well": well,
                        "File_Name": path.name,
                        "Status": "failed",
                        "Error": message,
                    }
                )
                records[well] = {
                    "status": "failed",
                    "error": message,
                    "stage": "Stitching"
                    if stitched_record["status"] != "completed"
                    else "Cell analysis",
                }
                log.record_error(f"Skipping {well}", error)
            checkpoint()
            log.phase(
                2,
                3,
                "Measure",
                finished=index + 1,
                count=len(metadata["wells"]),
            )
    except (Exception, KeyboardInterrupt) as error:
        failures += 1
        status["error"] = str(error) or type(error).__name__
        status["status"] = (
            "CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else ("NO_INPUT" if isinstance(error, NoInputError) else "FAILED")
        )
        log.record_error("Cell analysis", error)
        if isinstance(error, KeyboardInterrupt):
            raise
    finally:
        for well, record in records.items():
            if record["status"] == "pending":
                error = status.get("error", "folder interrupted")
                reason = f"Not completed: {error}"
                record.update(
                    status="failed", error=reason, stage="Cell analysis"
                )
                rows.append(
                    {
                        "Well": well,
                        "File_Name": f"{well}_stitched.tif",
                        "Status": "failed",
                        "Error": reason,
                    }
                )
        completed = sum(row["Status"] == "completed" for row in rows)
        failures = sum(
            row["status"] != "completed" for row in records.values()
        ) + bool(status.get("error"))
        status.update(
            status=status["status"]
            if status["status"] in {"CANCELLED", "NO_INPUT"}
            else outcome(completed, failures),
            completed_wells=completed,
            failures=failures,
            finished_utc=utc_now(),
        )
        try:
            log.phase(3, 3, "Save audit")
            save_exclusions(output, exclusions(status))
            checkpoint()
            log.event(
                status["status"],
                "Cell analysis",
                (
                    f"{completed} well(s) saved; {failures} failure(s). "
                    f"Results: {output}"
                ),
            )
        finally:
            log.close()
    return completed, failures


def run_analysis(args: argparse.Namespace) -> int:
    """Initialize one Fiji context for the batch and always dispose it."""
    folders = read_config(Path(args.input))
    context = None

    def ensure_context() -> None:
        nonlocal context
        if context is None:
            from scyjava import config, jvm_started

            if not jvm_started():
                config.add_option("-Xmx16g")
            context = initialize_imagej()

    completed, failures = 0, 0
    try:
        for folder in dict.fromkeys(folders):
            try:
                success, failed = process_folder(folder, args, ensure_context)
            except Exception as error:
                command_error("cell_count", error, folder)
                failures += 1
                continue
            completed += success
            failures += failed
    finally:
        if context is not None:
            context.dispose()
    print(
        f"Cell analysis: {completed} well(s) saved; {failures} failure(s).",
        flush=True,
    )
    return 1 if failures or not completed else 0


def main(argv: Sequence[str] | None = None) -> int:
    """Keep all Fiji worker cleanup at the standalone command boundary."""
    args = parse_args(argv)
    try:
        return _run_imagej_command(run_analysis, args)
    except KeyboardInterrupt as error:
        batch_error("cell_count", error, args.input)
        return 130
    except Exception as error:
        batch_error("cell_count", error, args.input)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
