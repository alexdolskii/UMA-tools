"""Segment prepared nuclei stacks without changing their preprocessing."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from importlib.metadata import version
from pathlib import Path

from uma_tools.config import folder_entries, load_json, resolve_folders
from uma_tools.files import save_csv, save_json, sha256_file
from uma_tools.run import scoped_file_log, utc_now

from .segment_models import bundled_models, choose_model

OUTPUT_DIRS = ("masks", "analysis", "clustering")


def read_settings(path: Path) -> tuple[dict, dict[Path, list[Path]]]:
    """Validate every source and output before loading the model."""
    raw = load_json(path)
    folders = list(
        dict.fromkeys(resolve_folders(folder_entries(raw), use_realpath=True))
    )
    tiles = raw.get("n_tiles", [1, 1, 1])
    if (
        not isinstance(tiles, (list, tuple))
        or len(tiles) != 3
        or any(type(value) is not int or value < 1 for value in tiles)
    ):
        raise ValueError("n_tiles must contain three positive integers: Z,Y,X")
    raw["n_tiles"] = tuple(tiles)
    inputs = {}
    for folder in folders:
        processed = folder / "processed"
        if not processed.is_dir():
            raise ValueError(
                f"Missing processed folder; run preparation: {processed}"
            )
        files = sorted(
            p
            for p in processed.iterdir()
            if p.is_file()
            and not p.name.startswith((".", "_"))
            and p.suffix.lower() in (".tif", ".tiff")
        )
        if not files:
            raise ValueError(f"No prepared TIFF stacks found: {processed}")
        names = [p.stem.casefold() for p in files]
        if len(set(names)) != len(names):
            raise ValueError(
                f"Prepared names would overwrite the same mask: {processed}"
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
                for source in folders
            ) or any(
                p.resolve().is_relative_to(output.resolve()) for p in files
            ):
                raise ValueError(f"An output contains an input: {output}")
        inputs[folder] = files
    all_files = [path.resolve() for files in inputs.values() for path in files]
    for folder in inputs:
        for name in OUTPUT_DIRS:
            output = (folder / name).resolve()
            if any(path.is_relative_to(output) for path in all_files):
                raise ValueError(f"An output contains an input: {output}")
    return raw, inputs


def confirm_replacement(inputs, overwrite: bool) -> None:
    """Avoid mixing new masks with stale masks or downstream calculations."""
    previous = [
        folder / name
        for folder in inputs
        for name in OUTPUT_DIRS
        if (folder / name).is_dir() and any((folder / name).iterdir())
    ]
    if not previous or overwrite:
        return
    print("Existing masks and dependent calculations will be replaced:")
    for path in previous:
        print(f"  {path}")
    if not sys.stdin.isatty():
        raise ValueError(
            "Previous results exist; use --overwrite to replace them"
        )
    if input("Replace these results? [y/N]: ").strip().lower() not in (
        "y",
        "yes",
    ):
        raise ValueError("Segmentation cancelled; previous results preserved")


def read_volume(path: Path):
    """Read one scalar ZYX stack and its saved calibration."""
    import numpy as np
    import tifffile

    with tifffile.TiffFile(path) as tif:
        if len(tif.series) != 1:
            raise ValueError(f"Expected one TIFF series: {path.name}")
        series = tif.series[0]
        array = series.asarray()
        axes = series.axes
        if array.ndim == 2 and axes == "YX":
            array = array[np.newaxis]
        elif array.ndim != 3 or axes not in ("ZYX", "QYX", "IYX"):
            raise ValueError(
                f"Expected scalar ZYX data, got {axes}: {path.name}"
            )
        if (
            array.dtype.kind not in "uif"
            or 0 in array.shape
            or not np.isfinite(array).all()
        ):
            raise ValueError(f"Invalid/nonfinite intensity data: {path.name}")
        metadata = dict(tif.imagej_metadata or {})
        calibration = {
            key: metadata[key]
            for key in ("spacing", "unit")
            if key in metadata
        }
        resolution = []
        for name in ("XResolution", "YResolution"):
            tag = tif.pages[0].tags.get(name)
            if tag is not None:
                numerator, denominator = tag.value
                if numerator > 0 and denominator > 0:
                    resolution.append((numerator, denominator))
        if len(resolution) == 2:
            calibration["resolution"] = resolution
        unit = tif.pages[0].tags.get("ResolutionUnit")
        if unit is not None:
            calibration["resolutionunit"] = int(unit.value)
    return array, calibration


def compact_dtype(labels):
    """Keep every object ID; never wrap labels greater than 65535."""
    import numpy as np

    if labels.dtype.kind not in "ui" or labels.min() < 0:
        raise ValueError("StarDist returned invalid instance labels")
    maximum = int(labels.max())
    if maximum > np.iinfo(np.uint32).max:
        raise ValueError("Instance labels exceed the supported uint32 range")
    return labels.astype(np.uint16 if maximum <= 65535 else np.uint32)


def predict(volume, model, tiles, logger):
    """Keep the original whole-volume 1/99.8 normalization and prediction."""
    import numpy as np
    from csbdeep.utils import normalize

    low, high = np.percentile(volume, (1, 99.8))
    if volume.min() == volume.max():
        logger.warning("Constant-intensity stack: writing an empty mask")
        return np.zeros(volume.shape, np.uint16), np.zeros(
            volume.shape, np.float32
        )
    if high <= low:
        raise ValueError(
            "1/99.8 percentile range is zero; inspect the sparse signal"
        )
    normalized = normalize(volume, 1, 99.8, axis=None)
    labels, _ = model.predict_instances(
        normalized,
        n_tiles=tiles,
        show_tile_progress=False,
    )
    if labels.shape != volume.shape:
        raise ValueError("Prediction changed the input stack dimensions")
    return compact_dtype(labels), normalized


def save_mask(path: Path, labels, calibration: dict) -> None:
    """Publish a complete calibrated TIFF, retaining large integer labels."""
    import tifffile

    temporary = path.with_suffix(".tmp.tif")
    metadata = {
        key: calibration[key]
        for key in ("spacing", "unit")
        if key in calibration
    }
    metadata["axes"] = "ZYX"
    options = {
        key: calibration[key]
        for key in ("resolution", "resolutionunit")
        if key in calibration
    }
    try:
        tifffile.imwrite(
            temporary,
            labels,
            imagej=labels.dtype.itemsize <= 2 and labels.shape[0] > 1,
            photometric="minisblack",
            metadata=metadata,
            **options,
        )
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def save_qc(path: Path, normalized, labels, count: int) -> None:
    """Show corresponding XY/XZ/YZ slices with stable object colors."""
    import matplotlib
    import numpy as np

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    z, y, x = (size // 2 for size in labels.shape)
    planes = [
        (normalized[z], labels[z]),
        (normalized[:, y], labels[:, y]),
        (normalized[:, :, x], labels[:, :, x]),
    ]
    ids = np.concatenate(([0], np.unique(labels[labels > 0])))
    colors = np.random.default_rng(6).uniform(0.2, 1.0, (len(ids), 4))
    colors[:, 3] = 0.45
    colors[0] = 0
    cmap = ListedColormap(colors)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
    try:
        for index, ((image, mask), name) in enumerate(
            zip(planes, ("XY", "XZ", "YZ"))
        ):
            for row in range(2):
                axes[row, index].imshow(image, cmap="gray", vmin=0, vmax=1)
                axes[row, index].set_title(
                    name if row == 0 else f"{name}: instances"
                )
                axes[row, index].axis("off")
            axes[1, index].imshow(
                np.searchsorted(ids, mask),
                cmap=cmap,
                vmin=0,
                vmax=max(len(ids) - 1, 1),
                interpolation="nearest",
            )
        fig.suptitle(
            f"{path.stem}: {count} nuclei; central slices (pixel axes)"
        )
        fig.savefig(path, dpi=120)
    finally:
        plt.close(fig)


def load_model(path: Path, thresholds: dict):
    """Load weights once, without initializing ImageJ or downloading models."""
    try:
        import numpy as np
        import tensorflow as tf
        from stardist.models import StarDist3D
    except ImportError as error:
        raise RuntimeError(
            "Install segmentation support in your UMA Python 3.10 "
            "environment: "
            'python -m pip install "./nuclei_layers_assay[segment]". '
            f"Dependency error: {error}"
        ) from error
    np.random.seed(6)
    model = StarDist3D(None, name=path.name, basedir=str(path.parent))
    if dict(model.thresholds._asdict()) != thresholds:
        raise ValueError("Loaded thresholds do not match thresholds.json")
    return model, tf


def run(inputs, settings: dict, model_path: Path, thresholds: dict) -> None:
    """Process all folders and expose every per-image failure to the caller."""
    import numpy as np

    model, tf = load_model(model_path, thresholds)
    logger = logging.Logger("uma_nla_segment", logging.INFO)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(console)
    total_failed = total_saved = 0
    try:
        provenance = {
            "model_path": str(model_path),
            "model_files_sha256": {
                name: sha256_file(model_path / name)
                for name in (
                    "config.json",
                    "thresholds.json",
                    "weights_best.h5",
                )
            },
            "thresholds": thresholds,
            "n_tiles": list(settings["n_tiles"]),
            "normalization": {
                "pmin": 1,
                "pmax": 99.8,
                "scope": "whole ZYX volume",
            },
            "versions": {
                name: version(name)
                for name in (
                    "uma-nuclei-layers-assay",
                    "uma-tools",
                    "stardist",
                    "csbdeep",
                    "tensorflow",
                    "numpy",
                )
            },
        }
        for folder, files in inputs.items():
            with scoped_file_log(logger, folder, "nuclei_analysis2.log"):
                logger.info("Source: %s; model: %s", folder, model_path)
                logger.info(
                    "Thresholds: %s; n_tiles: %s",
                    thresholds,
                    settings["n_tiles"],
                )
                if settings.get("downscale_factor") not in (None, 1):
                    logger.warning(
                        "Legacy downscale_factor is unused; "
                        "native resolution is preserved"
                    )
                for name in OUTPUT_DIRS:
                    output = folder / name
                    if output.exists():
                        shutil.rmtree(output)
                    output.mkdir()
                masks = folder / "masks"
                status = {
                    "status": "started",
                    "started_utc": utc_now(),
                    **provenance,
                    "images": [],
                }
                save_json(masks / "run_status.json", status)
                completed = False
                try:
                    for source in files:
                        mask_path = masks / f"{source.stem}_mask.tif"
                        qc_path = masks / f"{source.stem}_QC.png"
                        record = {"file": source.name, "status": "failed"}
                        try:
                            volume, calibration = read_volume(source)
                            labels, normalized = predict(
                                volume, model, settings["n_tiles"], logger
                            )
                            count = int(np.count_nonzero(np.unique(labels)))
                            save_mask(mask_path, labels, calibration)
                            save_qc(qc_path, normalized, labels, count)
                            record.update(
                                status="success",
                                nuclei=count,
                                shape=list(labels.shape),
                                dtype=str(labels.dtype),
                                calibration=calibration,
                                input_sha256=sha256_file(source),
                            )
                        except Exception as error:
                            total_failed += 1
                            mask_path.unlink(missing_ok=True)
                            qc_path.unlink(missing_ok=True)
                            record["error"] = str(error)
                            logger.exception("FAILED: %s", source.name)
                        else:
                            total_saved += 1
                            logger.info(
                                "Saved %s: %s nuclei; %s",
                                mask_path.name,
                                count,
                                labels.shape,
                            )
                        status["images"].append(record)
                    completed = True
                finally:
                    status["status"] = (
                        "success"
                        if completed
                        and all(
                            record["status"] == "success"
                            for record in status["images"]
                        )
                        else "failed"
                    )
                    status["finished_utc"] = utc_now()
                    save_json(masks / "run_status.json", status)
                    save_csv(
                        masks / "segmentation_summary.csv",
                        ["File_Name", "Status", "Nuclei", "Error"],
                        [
                            {
                                "File_Name": r["file"],
                                "Status": r["status"],
                                "Nuclei": r.get("nuclei", ""),
                                "Error": r.get("error", ""),
                            }
                            for r in status["images"]
                        ],
                    )
    finally:
        logger.removeHandler(console)
        console.close()
        tf.keras.backend.clear_session()
    print(
        f"Nuclei segmentation: {total_saved} saved; {total_failed} failed.",
        flush=True,
    )
    if total_failed:
        raise RuntimeError("Some images failed; see nuclei_analysis2.log")


def main(argv=None, *, default_input: str | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="StarDist 3D nuclei segmentation"
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
        default=default_input,
        help="The same nuclei_layers.json used for preparation",
    )
    parser.add_argument(
        "--model", help="Bundled model name or custom model directory"
    )
    parser.add_argument(
        "--list-models",
        action="store_true",
        help="List bundled models without loading TensorFlow",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace masks, analysis and clustering without a prompt",
    )
    args = parser.parse_args(argv)
    if args.list_models:
        for name in bundled_models():
            print(name)
        return 0
    if args.input is None:
        parser.error("-i/--input is required for segmentation")
    try:
        settings, inputs = read_settings(args.input)
        model_path, thresholds = choose_model(
            args.model
            if args.model is not None
            else settings.get("model_path"),
            args.input,
            sys.stdin.isatty(),
        )
        for folder in inputs:
            if any(
                model_path.is_relative_to(folder / name)
                for name in OUTPUT_DIRS
            ):
                raise ValueError("The model cannot be inside an output folder")
        confirm_replacement(inputs, args.overwrite)
        run(inputs, settings, model_path, thresholds)
        return 0
    except KeyboardInterrupt:
        print("Nuclei segmentation interrupted.", file=sys.stderr, flush=True)
        return 130
    except Exception as error:
        print(
            f"Nuclei segmentation failed: {error}", file=sys.stderr, flush=True
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
