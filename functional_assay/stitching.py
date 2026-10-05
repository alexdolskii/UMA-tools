"""Stitch nine single-channel ND2 stacks per well for cell counting.

Preserve the supplied grid, fusion parameters, full stacks, Sharpen
filter, and TIFF outputs. The installed command owns worker shutdown.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections import Counter
from collections.abc import Sequence
from importlib.metadata import version
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from uma_tools.cli import _run_imagej_command
from uma_tools.config import read_config
from uma_tools.files import save_json, sha256_file
from uma_tools.imagej import FIJI_ENDPOINT, initialize_imagej
from uma_tools.run import utc_now

from .calibration import read_well_calibration
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
    warning,
)

FRAME_POSITIONS = {0: 1, 1: 2, 2: 3, 5: 4, 4: 5, 3: 6, 6: 7, 7: 8, 8: 9}
FILENAME_PATTERN = re.compile(
    r".*__(Well[A-Za-z0-9]+)_Point[A-Za-z0-9]+_(\d{4})_Channel.*\.nd2$",
    re.IGNORECASE,
)


def discover_wells(folder: Path) -> dict[str, list[tuple[int, Path]]]:
    """Group visible ND2 files using the original well/frame pattern."""
    wells: dict[str, list[tuple[int, Path]]] = {}
    for path in sorted(folder.iterdir()):
        if path.name.startswith(".") or not path.is_file():
            continue
        if path.suffix.lower() != ".nd2":
            continue
        match = FILENAME_PATTERN.fullmatch(path.name)
        if match is None:
            warning(f"Skipping unrecognized ND2 filename: {path.name}")
            continue
        well_id, frame = match.groups()
        wells.setdefault(well_id, []).append((int(frame), path))
    return wells


def validate_frames(files: list[tuple[int, Path]]) -> None:
    """Require exactly one file for each of the nine expected frames."""
    counts = Counter(index for index, _ in files)
    expected = set(FRAME_POSITIONS)
    missing = sorted(expected - counts.keys())
    unexpected = sorted(counts.keys() - expected)
    duplicates = sorted(index for index, count in counts.items() if count > 1)
    if len(files) != 9 or missing or unexpected or duplicates:
        raise ValueError(
            "Expected nine unique frames 0000-0008; "
            f"found {len(files)} files, missing={missing}, "
            f"unexpected={unexpected}, duplicates={duplicates}"
        )


def build_stitching_options(folder: Path, overlap: float) -> str:
    """Keep the supplied plugin settings and parameterize overlap."""
    directory = folder.resolve().as_posix() + "/"
    if "]" in directory:
        raise ValueError("ImageJ stitching paths cannot contain ']'")
    return (
        "type=[Grid: row-by-row] "
        "order=[Right & Down] "
        "grid_size_x=3 grid_size_y=3 "
        f"tile_overlap={overlap} "
        "first_file_index_i=1 "
        f"directory=[{directory}] "
        "file_names=[image_{i}.nd2] "
        "output_textfile_name=TileConfiguration.txt "
        "fusion_method=[Linear Blending] "
        "regression_threshold=0.30 "
        "max/avg_displacement_threshold=2.50 "
        "absolute_displacement_threshold=3.50 "
        "computation_parameters=[Save computation time (but use more RAM)] "
        "image_output=[Fuse and display]"
    )


def stitch_well(
    folder: Path,
    well_id: str,
    files: list[tuple[int, Path]],
    output_folder: Path,
    overlap: float,
    record: dict | None = None,
) -> Path:
    """Fuse full stacks in batch mode, sharpen, and save one TIFF."""
    from scyjava import jimport

    validate_frames(files)
    ij = jimport("ij.IJ")
    interpreter = jimport("ij.macro.Interpreter")
    with TemporaryDirectory(
        prefix=f".uma_stitch_{well_id}_", dir=folder
    ) as tmp:
        temporary = Path(tmp)
        for index, source in files:
            destination = temporary / f"image_{FRAME_POSITIONS[index]}.nd2"
            try:
                destination.symlink_to(source.resolve())
            except OSError:
                shutil.copy2(source, destination)

        options = build_stitching_options(temporary, overlap)
        quoted_options = json.dumps(options, ensure_ascii=False)
        macro = f'run("Grid/Collection stitching", {quoted_options});'
        image = interpreter().runBatchMacro(macro, None)
        if image is None:
            raise RuntimeError(f"No stitched image was produced for {well_id}")
        try:
            if image.getNChannels() != 1:
                raise ValueError(f"{well_id}: exactly one channel is required")
            activity(f"{well_id}: Sharpen, all slices")
            ij.run(image, "Sharpen", "stack")
            output_file = output_folder / f"{well_id}_stitched.tif"
            ij.saveAs(image, "Tiff", str(output_file))
            if not output_file.is_file() or output_file.stat().st_size == 0:
                raise OSError(f"Could not save stitched TIFF: {output_file}")
            if record is not None:
                record.update(
                    output_file=output_file.name,
                    width_px=int(image.getWidth()),
                    height_px=int(image.getHeight()),
                    slices=int(image.getNSlices()),
                    sha256=sha256_file(output_file),
                )
            activity(f"{well_id}: TIFF saved (overlap {overlap}%)")
            return output_file
        finally:
            image.close()


def reset_output_folder(folder: Path) -> Path:
    """Replace the previous Stitched_Results directory as requested."""
    output = assay_directory(folder) / "Stitched_Results"
    if output.is_symlink():
        raise ValueError(
            f"Output directory must not be a symbolic link: {output}"
        )
    if output.exists():
        shutil.rmtree(output)
    output.mkdir()
    return output


def process_folder(folder: Path, overlap: float, ensure_context=None) -> dict:
    """Replace stale outputs and finish every independent well."""
    output = reset_output_folder(folder)
    log = RunLog(output, folder, "stitching")
    metadata = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "functional_assay_version": version("uma-functional-assay"),
        "uma_tools_version": version("uma-tools"),
        "fiji_endpoint": FIJI_ENDPOINT,
        "overlap_percent": overlap,
        "grid": "3x3, row-by-row, right and down",
        "status": "RUNNING",
        "source": str(folder),
        "output": str(output),
        "wells": {},
    }
    metadata_path = output / "stitching_metadata.json"
    try:
        save_json(metadata_path, metadata, allow_nan=False)
        log.event("STARTED", "Stitching", folder)
        log.event("INFO", "Output", output)
        log.event(
            "INFO", "Plan", "1/3 discover; 2/3 stitch wells; 3/3 save audit"
        )
        log.phase(1, 3, "Discover frames")
        wells = discover_wells(folder)
        if not wells:
            raise NoInputError("No eligible ND2 wells found")
        metadata["wells"] = {well: {"status": "pending"} for well in wells}
        for index, (well_id, files) in enumerate(wells.items()):
            log.phase(
                2,
                3,
                "Stitch",
                finished=index,
                count=len(wells),
                detail=f"{well_id}: validate nine frames",
            )
            record = {
                "status": "running",
                "source_frames": [
                    {
                        "filename": path.name,
                        "frame_index": index,
                        "grid_position": FRAME_POSITIONS.get(index),
                    }
                    for index, path in sorted(files)
                ],
            }
            metadata["wells"][well_id] = record
            try:
                validate_frames(files)
                if ensure_context is not None:
                    activity(f"{well_id}: initialize ImageJ if needed")
                    ensure_context()
                try:
                    record["calibration"] = read_well_calibration(files)
                except Exception as error:
                    # The TIFF remains useful; cell analysis requires a scale.
                    record["calibration_error"] = str(error)
                    log.event("WARNING", well_id, f"Calibration: {error}")
                activity(f"{well_id}: fuse nine stacks")
                stitch_well(folder, well_id, files, output, overlap, record)
                record["status"] = "completed"
                log.event(
                    "INFO", well_id, "Stitched TIFF saved", console=False
                )
            except Exception as error:
                record.update(
                    status="failed",
                    error=str(error) or type(error).__name__,
                    stage="Stitching",
                )
                try:
                    (output / f"{well_id}_stitched.tif").unlink(
                        missing_ok=True
                    )
                except OSError as cleanup_error:
                    log.event(
                        "WARNING",
                        well_id,
                        "Could not remove incomplete TIFF: "
                        + str(cleanup_error),
                    )
                log.record_error(f"Skipping {well_id}", error)
            finally:
                save_json(metadata_path, metadata, allow_nan=False)
            log.phase(2, 3, "Stitch", finished=index + 1, count=len(wells))
        completed = sum(
            r["status"] == "completed" for r in metadata["wells"].values()
        )
        metadata["status"] = outcome(completed, len(wells) - completed)
    except (Exception, KeyboardInterrupt) as error:
        metadata.update(
            status="CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else ("NO_INPUT" if isinstance(error, NoInputError) else "FAILED"),
            error=str(error) or type(error).__name__,
        )
        log.record_error("Stitching", error)
        if isinstance(error, KeyboardInterrupt):
            raise
    finally:
        records = metadata["wells"]
        for record in records.values():
            if record["status"] in {"pending", "running"}:
                error = metadata.get("error", metadata["status"])
                record.update(
                    status="failed",
                    stage="Stitching",
                    error=f"Not completed: {error}",
                )
        metadata["completed_wells"] = sum(
            r["status"] == "completed" for r in records.values()
        )
        metadata["failures"] = len(records) - metadata["completed_wells"]
        metadata["finished_utc"] = utc_now()
        try:
            log.phase(3, 3, "Save audit")
            save_exclusions(output, exclusions(metadata))
            save_json(metadata_path, metadata, allow_nan=False)
            save_json(output / "run_status.json", metadata, allow_nan=False)
            log.event(
                metadata["status"],
                "Stitching",
                (
                    f"{metadata['completed_wells']} well(s) saved; "
                    f"{metadata['failures']} excluded. Results: {output}"
                ),
            )
        finally:
            log.close()
    return metadata


def process_wells_stitching(json_path: str, overlap: float = 32.8) -> int:
    """Read UMA input folders, initialize Fiji once, and dispose it."""
    if not math.isfinite(overlap) or not 0 <= overlap < 100:
        raise ValueError("Overlap must be a finite percentage from 0 to <100")
    folders = read_config(Path(json_path))

    from scyjava import config, jvm_started

    # Preserve the original 16 GiB heap cap without global Java flags.
    if not jvm_started():
        config.add_option("-Xmx16g")
    context: Any = None

    def ensure_context():
        nonlocal context
        if context is None:
            context = initialize_imagej()

    results = []
    try:
        for folder in dict.fromkeys(folders):
            try:
                results.append(process_folder(folder, overlap, ensure_context))
            except Exception as error:
                command_error("stitching", error, folder)
                results.append({"status": "ERROR"})
    finally:
        if context is not None:
            context.dispose()
    return (
        0 if results and all(r["status"] == "SUCCESS" for r in results) else 1
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the optional stitching command in the existing UMA environment."""
    parser = argparse.ArgumentParser(
        description="Stitch nine ND2 stacks per functional-assay well"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('uma-functional-assay')}",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to a JSON file containing folder_paths",
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=32.8,
        help="Tile overlap in percent, from 0 to below 100 (default: 32.8)",
    )
    args = parser.parse_args(argv)
    if not math.isfinite(args.overlap) or not 0 <= args.overlap < 100:
        parser.error("--overlap must be a finite percentage from 0 to <100")
    try:
        return _run_imagej_command(
            process_wells_stitching, args.input, args.overlap
        )
    except KeyboardInterrupt as error:
        batch_error("stitching", error, args.input)
        return 130
    except Exception as error:
        batch_error("stitching", error, args.input)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
