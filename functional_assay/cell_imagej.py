"""Fiji preprocessing, area masks, Watershed counts, and contour images."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage


def filter_particles(
    mask: np.ndarray, minimum_px: float, *, exclude_edges: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """Keep 8-connected objects; preserve holes and optionally omit edges.

    Return consecutive integer labels and their foreground pixel counts.
    No hole filling is performed, so areas measure positive mask pixels.
    """
    labels, count = ndimage.label(mask, structure=np.ones((3, 3)))
    sizes = np.bincount(labels.ravel(), minlength=count + 1)
    keep = sizes >= minimum_px
    keep[0] = False
    if exclude_edges:
        edges = np.concatenate(
            (labels[0], labels[-1], labels[:, 0], labels[:, -1])
        )
        keep[np.unique(edges)] = False
    indices = np.flatnonzero(keep)
    mapping = np.zeros(len(sizes), dtype=np.int32)
    mapping[indices] = np.arange(1, len(indices) + 1)
    return mapping[labels], sizes[indices]


def mask_image(mask: np.ndarray, title: str) -> Any:
    """Create an ImageJ binary image, with foreground 255 and background 0."""
    import jpype
    from scyjava import jimport

    height, width = mask.shape
    pixels = (mask.astype(np.uint8) * 255).view(np.int8).ravel()
    processor = jimport("ij.process.ByteProcessor")(
        width, height, jpype.JArray(jpype.JByte)(pixels), None
    )
    return jimport("ij.ImagePlus")(title, processor)


def save_mask(path: Path, mask: np.ndarray, calibration: dict) -> None:
    """Write a calibrated, binary TIFF without altering the input image."""
    from scyjava import jimport

    image = mask_image(mask, path.stem)
    try:
        scale = jimport("ij.measure.Calibration")()
        scale.pixelWidth = calibration["pixel_size_x_um"]
        scale.pixelHeight = calibration["pixel_size_y_um"]
        scale.setUnit("um")
        image.setCalibration(scale)
        if not jimport("ij.io.FileSaver")(image).saveAsTiff(str(path)):
            raise OSError(f"Could not save mask: {path}")
    finally:
        image.close()


def prepare_projection(path: Path) -> Any:
    """Apply the supplied macros' MAX, background, contrast, and smoothing."""
    from scyjava import jimport

    ij = jimport("ij.IJ")
    source = ij.openImage(str(path.resolve()))
    if source is None:
        raise ValueError(f"Could not open stitched TIFF: {path}")
    projection = None
    try:
        if source.getNChannels() != 1 or source.getNFrames() != 1:
            raise ValueError("Stitched TIFF must have one channel/time point")
        if source.getBitDepth() != 16:
            raise ValueError("Stitched TIFF must contain 16-bit intensities")
        if source.getStackSize() > 1:
            projection = jimport("ij.plugin.ZProjector").run(source, "max")
        else:
            projection = source.duplicate()
        # All algorithm sizes are pixels. Physical scale is applied to output.
        projection.setCalibration(jimport("ij.measure.Calibration")())
        ij.run(projection, "Subtract Background...", "rolling=50")
        ij.run(projection, "Enhance Contrast", "saturated=0.35")
        ij.run(projection, "Smooth", "")
        ij.run(projection, "Smooth", "")
        return projection
    except BaseException:
        if projection is not None:
            projection.close()
        raise
    finally:
        source.close()


def analyze_image(
    path: Path, thresholds: tuple[float, float] | None, minimum_px: float
) -> dict[str, Any]:
    """Measure cleaned area before Watershed; count accepted objects after."""
    from scyjava import jimport

    ij = jimport("ij.IJ")
    projection = prepare_projection(path)
    try:
        processor = projection.getProcessor()
        if thresholds is None:
            ij.setAutoThreshold(projection, "RenyiEntropy dark no-reset")
        else:
            ij.setThreshold(projection, *thresholds)
        lower = float(processor.getMinThreshold())
        upper = float(processor.getMaxThreshold())
        if not 0 <= lower <= upper <= 65535:
            raise ValueError(f"Invalid effective threshold: {lower}, {upper}")
        height, width = projection.getHeight(), projection.getWidth()
        raw_mask = (
            np.asarray(processor.createMask().getPixels()).reshape(
                height, width
            )
            != 0
        )
        labels, _ = filter_particles(raw_mask, minimum_px)
        area_mask = labels != 0
        separated = mask_image(area_mask, "Watershed")
        preferences = jimport("ij.Prefs")
        old_background = preferences.blackBackground
        try:
            preferences.blackBackground = True
            ij.run(separated, "Watershed", "")
            divided = (
                np.asarray(separated.getProcessor().getPixels()).reshape(
                    height, width
                )
                != 0
            )
            object_labels, areas = filter_particles(
                divided, minimum_px, exclude_edges=True
            )
        finally:
            preferences.blackBackground = old_background
            separated.close()
        values = np.asarray(processor.getPixels()).astype(np.uint16)
        values = values.reshape(height, width)
        return {
            "area_mask": area_mask,
            "counting_mask": object_labels != 0,
            "object_areas_px2": areas,
            "projection": values,
            "display_min": float(projection.getDisplayRangeMin()),
            "display_max": float(projection.getDisplayRangeMax()),
            "threshold_lower": lower,
            "threshold_upper": upper,
            "width_px": int(width),
            "height_px": int(height),
        }
    finally:
        projection.close()


def save_contours(path: Path, result: dict[str, Any]) -> None:
    """Overlay counted-object contours in green on the processed MAX image."""
    from PIL import Image

    values = result["projection"].astype(float)
    low, high = result["display_min"], result["display_max"]
    scaled = np.clip((values - low) / max(high - low, 1), 0, 1)
    gray = (scaled * 255).astype(np.uint8)
    rgb = np.repeat(gray[:, :, None], 3, axis=2)
    counted = result["counting_mask"]
    contours = counted & ~ndimage.binary_erosion(
        counted, structure=np.ones((3, 3)), border_value=0
    )
    rgb[contours] = (0, 255, 0)
    Image.fromarray(rgb).save(path)
