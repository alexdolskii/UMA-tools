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
from uma_tools.imagej import initialize_imagej

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
            print(f"Skipping unrecognized ND2 filename: {path.name}")
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
            print("  Applying Sharpen to all slices in the stitched stack...")
            ij.run(image, "Sharpen", "stack")
            output_file = output_folder / f"{well_id}_stitched.tif"
            ij.saveAs(image, "Tiff", str(output_file))
            if not output_file.is_file() or output_file.stat().st_size == 0:
                raise OSError(f"Could not save stitched TIFF: {output_file}")
            print(f"Saved: {output_file} (overlap: {overlap}%)")
            return output_file
        finally:
            image.close()


def reset_output_folder(folder: Path) -> Path:
    """Replace the previous Stitched_Results directory as requested."""
    output = folder / "Stitched_Results"
    if output.is_symlink():
        raise ValueError(
            f"Output directory must not be a symbolic link: {output}"
        )
    if output.exists():
        print(f"Replacing previous results: {output}")
        shutil.rmtree(output)
    output.mkdir()
    return output


def process_folder(folder: Path, overlap: float) -> int:
    """Skip invalid frame sets and stitch the complete wells."""
    if not folder.is_dir():
        print(f"Folder not found, skipping: {folder}")
        return 0
    print(f"\n--- Processing folder: {folder} ---")
    wells = discover_wells(folder)
    valid_wells = {}
    for well_id, files in wells.items():
        try:
            validate_frames(files)
        except ValueError as error:
            print(f"Skipping {well_id}: {error}")
            continue
        valid_wells[well_id] = files
    if not valid_wells:
        print("No complete nine-frame wells found, skipping folder.")
        return 0

    output = reset_output_folder(folder)
    for well_id, files in valid_wells.items():
        print(f"Stitching {well_id}...")
        stitch_well(folder, well_id, files, output, overlap)
    return len(valid_wells)


def process_wells_stitching(json_path: str, overlap: float = 32.8) -> None:
    """Read UMA input folders, initialize Fiji once, and dispose it."""
    if not math.isfinite(overlap) or not 0 <= overlap < 100:
        raise ValueError("Overlap must be a finite percentage from 0 to <100")
    folders = read_config(Path(json_path))

    from scyjava import config, jvm_started

    # Preserve the original 16 GiB heap cap without global Java flags.
    if not jvm_started():
        config.add_option("-Xmx16g")
    context: Any = initialize_imagej()
    try:
        count = sum(process_folder(folder, overlap) for folder in folders)
        if count == 0:
            raise ValueError("No complete wells were stitched")
        print(f"Stitching completed: {count} well(s).")
    finally:
        context.dispose()


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
    return _run_imagej_command(
        process_wells_stitching, args.input, args.overlap
    )


if __name__ == "__main__":
    raise SystemExit(main())
