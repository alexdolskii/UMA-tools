"""Extract one nuclei channel and apply the original 3D filters."""

from __future__ import annotations

import argparse
import logging
import math
import shutil
import sys
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

from uma_tools.cli import _run_imagej_command
from uma_tools.config import folder_entries, load_json, resolve_folders
from uma_tools.imagej import FIJI_ENDPOINT, initialize_imagej

OUTPUT_DIRS = ("processed", "masks", "analysis", "clustering")


@dataclass(frozen=True)
class Settings:
    """Preparation parameters; later-stage JSON keys are ignored."""

    folders: tuple[Path, ...]
    channel: int
    gaussian_sigma: float
    mean_radius: int


def read_settings(path: Path) -> Settings:
    """Validate sources before starting Fiji or removing outputs."""
    raw = load_json(path)
    folders = tuple(
        dict.fromkeys(resolve_folders(folder_entries(raw), use_realpath=True))
    )
    for folder in folders:
        if not folder.is_dir():
            raise ValueError(f"Source folder does not exist: {folder}")
    channel = raw.get("nuclei_channel", 1)
    radius = raw.get("mean_radius", 3)
    sigma = raw.get("gaussian_sigma", 4.0)
    if type(channel) is not int or channel < 1:
        raise ValueError("nuclei_channel must be an integer starting from 1")
    if type(radius) is not int or radius < 0:
        raise ValueError("mean_radius must be a non-negative integer")
    if (
        type(sigma) not in (int, float)
        or not math.isfinite(sigma)
        or sigma < 0
    ):
        raise ValueError("gaussian_sigma must be finite and non-negative")
    return Settings(folders, channel, float(sigma), radius)


def discover_inputs(settings: Settings) -> dict[Path, list[Path]]:
    """Ignore metadata and reject ambiguous or unsafe output paths."""
    inputs = {}
    for folder in settings.folders:
        files = sorted(
            path
            for path in folder.iterdir()
            if path.is_file()
            and not path.name.startswith((".", "_"))
            and path.suffix.lower() in (".nd2", ".tif", ".tiff")
        )
        if not files:
            raise ValueError(f"No ND2/TIFF source images found: {folder}")
        names = [path.stem.casefold() for path in files]
        if len(set(names)) != len(names):
            raise ValueError(
                f"Source names would overwrite the same nuclei TIFF: {folder}"
            )
        for name in OUTPUT_DIRS:
            output = folder / name
            if output.is_symlink() or (
                output.exists() and not output.is_dir()
            ):
                raise ValueError(
                    f"Output must be a regular directory: {output}"
                )
            if any(
                source.resolve().is_relative_to(output.resolve())
                for source in settings.folders
            ):
                raise ValueError(f"Output contains a source folder: {output}")
        inputs[folder] = files
    return inputs


def confirm_replacement(settings: Settings, overwrite: bool) -> None:
    """Keep the original reset policy, with one batch confirmation."""
    existing = [
        folder / name
        for folder in settings.folders
        for name in OUTPUT_DIRS
        if (folder / name).exists()
    ]
    if not existing or overwrite:
        return
    print("The following previous results will be removed:", flush=True)
    for path in existing:
        print(f"  {path}", flush=True)
    if not sys.stdin.isatty():
        raise ValueError(
            "Previous results exist. Use --overwrite to replace them."
        )
    reply = input("Replace these results? [y/N]: ").strip().lower()
    if reply not in ("y", "yes"):
        raise ValueError("Preparation cancelled; previous results preserved")


class NucleiProcessor:
    """ImageJ object operations without a current image or GUI windows."""

    def __init__(self) -> None:
        from scyjava import jimport

        self.IJ = jimport("ij.IJ")
        self.Duplicator = jimport("ij.plugin.Duplicator")
        self.FileSaver = jimport("ij.io.FileSaver")
        self.BF = jimport("loci.plugins.BF")
        self.Options = jimport("loci.plugins.in.ImporterOptions")
        self.ImageReader = jimport("loci.formats.ImageReader")

    @staticmethod
    def release(image) -> None:
        if image is not None:
            image.changes = False
            image.close()

    def open_image(self, path: Path, channel: int):
        """Import a complete ND2 field or an ImageJ TIFF hyperstack."""
        if path.suffix.lower() != ".nd2":
            image = self.IJ.openImage(str(path))
            if image is None:
                raise ValueError(f"ImageJ could not open TIFF: {path.name}")
            return image
        reader = self.ImageReader()
        try:
            reader.setGroupFiles(False)
            reader.setId(str(path))
            if reader.getSeriesCount() != 1:
                raise ValueError(
                    f"Export multi-series ND2 as separate fields: {path.name}"
                )
            if reader.getSizeT() != 1 or reader.isRGB():
                raise ValueError(
                    f"Expected one time point and scalar channels: {path.name}"
                )
            if channel > reader.getSizeC():
                raise ValueError(
                    f"Channel {channel} exceeds {reader.getSizeC()} "
                    f"channels in {path.name}"
                )
            expected = [
                reader.getSizeX(),
                reader.getSizeY(),
                reader.getSizeC(),
                reader.getSizeZ(),
                reader.getSizeT(),
            ]
        finally:
            reader.close()
        options = self.Options()
        options.setId(str(path))
        # Match Bio-Formats' original default before conversion to 8-bit.
        options.setAutoscale(True)
        options.setQuiet(True)
        options.setWindowless(True)
        options.setGroupFiles(False)
        options.setCrop(False)
        options.setSpecifyRanges(False)
        options.setSplitChannels(False)
        options.setSplitFocalPlanes(False)
        options.setSplitTimepoints(False)
        options.setSwapDimensions(False)
        options.setConcatenate(False)
        options.setVirtual(False)
        options.setShowMetadata(False)
        options.setShowOMEXML(False)
        options.setShowROIs(False)
        options.setColorMode(self.Options.COLOR_MODE_DEFAULT)
        options.setStackFormat(self.Options.VIEW_HYPERSTACK)
        options.setStackOrder(self.Options.ORDER_XYCZT)
        options.setOpenAllSeries(False)
        options.clearSeries()
        options.setSeriesOn(0, True)
        images = self.BF.openImagePlus(options)
        if images is None or len(images) != 1:
            for image in images or []:
                self.release(image)
            raise ValueError(f"Expected one imported field: {path.name}")
        image = images[0]
        if list(image.getDimensions()) != expected:
            self.release(image)
            raise ValueError(f"ND2 import changed the dimensions: {path.name}")
        return image

    def process(self, source: Path, target: Path, settings: Settings) -> str:
        """Preserve Z/calibration and the original 8-bit/filter sequence."""
        original = filtered = None
        try:
            original = self.open_image(source, settings.channel)
            width, height, channels, slices, frames = original.getDimensions()
            if (
                frames != 1
                or not 1 <= settings.channel <= channels
                or original.getBitDepth() not in (8, 16, 32)
                or original.getStackSize() != channels * slices
            ):
                raise ValueError(
                    f"Unsupported channel, time points or dimensions in "
                    f"{source.name}: {list(original.getDimensions())}"
                )
            filtered = self.Duplicator().run(
                original,
                settings.channel,
                settings.channel,
                1,
                slices,
                1,
                1,
            )
            filtered.setCalibration(original.getCalibration().copy())
            self.IJ.run(filtered, "8-bit", "")
            sigma = settings.gaussian_sigma
            self.IJ.run(
                filtered,
                "Gaussian Blur 3D...",
                f"x={sigma} y={sigma} z={sigma}",
            )
            radius = settings.mean_radius
            self.IJ.run(
                filtered, "Mean 3D...", f"x={radius} y={radius} z={radius}"
            )
            saver = self.FileSaver(filtered)
            saved = (
                saver.saveAsTiffStack(str(target))
                if filtered.getStackSize() > 1
                else saver.saveAsTiff(str(target))
            )
            if not saved:
                raise OSError(f"ImageJ could not save TIFF: {target}")
            calibration = filtered.getCalibration()
            return (
                f"{width} x {height} x {slices} Z; 8-bit; "
                f"pixel size {calibration.pixelWidth:g} x "
                f"{calibration.pixelHeight:g} x {calibration.pixelDepth:g} "
                f"{calibration.getUnit()}"
            )
        finally:
            self.release(filtered)
            self.release(original)


def run(settings: Settings, inputs: dict[Path, list[Path]]) -> None:
    """Use one Fiji context; continue after an individual image failure."""
    from uma_tools.run import scoped_file_log

    logger = logging.getLogger("nuclei_layers_assay.prepare")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(console)
    context = None
    failed = saved = 0
    try:
        context = initialize_imagej()
        processor = NucleiProcessor()
        for folder, files in inputs.items():
            with scoped_file_log(logger, folder, "nuclei_analysis1.log"):
                logger.info("Source: %s; Fiji: %s", folder, FIJI_ENDPOINT)
                logger.info(
                    "Nuclei channel=%s; Gaussian sigma=%s; mean radius=%s",
                    settings.channel,
                    settings.gaussian_sigma,
                    settings.mean_radius,
                )
                for name in OUTPUT_DIRS:
                    output = folder / name
                    if output.exists():
                        shutil.rmtree(output)
                    output.mkdir()
                for source in files:
                    target = folder / "processed" / f"{source.stem}_nuclei.tif"
                    try:
                        details = processor.process(source, target, settings)
                    except Exception:
                        failed += 1
                        target.unlink(missing_ok=True)
                        logger.exception("FAILED: %s", source.name)
                    else:
                        saved += 1
                        logger.info("Saved %s; %s", target.name, details)
                logger.info("Folder processing finished: %s", folder)
    finally:
        try:
            if context is not None:
                active_error = sys.exc_info()[0] is not None
                try:
                    context.dispose()
                except Exception:
                    if not active_error:
                        raise
                    logger.exception("Could not dispose ImageJ after failure")
        finally:
            logger.removeHandler(console)
            console.close()
    print(f"Nuclei preparation: {saved} saved; {failed} failed.", flush=True)
    if failed:
        raise RuntimeError("Some images failed; see nuclei_analysis1.log")


def main(argv=None, *, default_input: str | None = None) -> int:
    """Parse options before importing or starting the scientific runtime."""
    parser = argparse.ArgumentParser(
        description="Prepare the nuclei channel with the original 3D filters"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {version('uma-nuclei-layers-assay')}",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=default_input is None,
        default=default_input,
        help="Path to nuclei_layers.json",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace processed, masks, analysis and clustering "
        "without a prompt",
    )
    args = parser.parse_args(argv)
    try:
        settings = read_settings(args.input)
        inputs = discover_inputs(settings)
        confirm_replacement(settings, args.overwrite)
        return _run_imagej_command(run, settings, inputs)
    except KeyboardInterrupt:
        print("Nuclei preparation interrupted.", file=sys.stderr, flush=True)
        return 130
    except Exception as error:
        print(
            f"Nuclei preparation failed: {error}", file=sys.stderr, flush=True
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
