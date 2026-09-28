"""Read and validate physical pixel sizes from nine original ND2 tiles."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any


def read_nd2_metadata(path: Path) -> dict[str, Any]:
    """Read OME metadata without loading image pixels or opening windows."""
    from scyjava import jimport

    service = jimport("loci.common.services.ServiceFactory")().getInstance(
        jimport("loci.formats.services.OMEXMLService")
    )
    metadata = service.createOMEXMLMetadata()
    reader = jimport("loci.formats.ImageReader")()
    reader.setMetadataStore(metadata)
    try:
        reader.setId(str(path.resolve()))
        if reader.getSeriesCount() != 1:
            raise ValueError(f"{path.name}: exactly one image series required")
        micrometer = jimport("ome.units.UNITS").MICROMETER
        sizes = {}
        for axis in ("X", "Y"):
            quantity = getattr(metadata, f"getPixelsPhysicalSize{axis}")(0)
            value = None if quantity is None else quantity.value(micrometer)
            if value is None:
                raise ValueError(f"{path.name}: missing physical size {axis}")
            size = float(value.doubleValue())
            if not math.isfinite(size) or size <= 0:
                raise ValueError(f"{path.name}: invalid physical size {axis}")
            sizes[f"pixel_size_{axis.lower()}_um"] = size
        return {
            "filename": path.name,
            "width_px": int(reader.getSizeX()),
            "height_px": int(reader.getSizeY()),
            "slices": int(reader.getSizeZ()),
            "channels": int(reader.getSizeC()),
            "frames": int(reader.getSizeT()),
            **sizes,
        }
    finally:
        reader.close()


def read_well_calibration(
    files: list[tuple[int, Path]],
) -> dict[str, Any]:
    """Require nine consistent, single-channel tiles; report sizes in µm."""
    indices = [index for index, _ in files]
    if sorted(indices) != list(range(9)):
        raise ValueError("Calibration requires nine unique frames 0000-0008")
    tiles = []
    for index, path in sorted(files):
        tile = read_nd2_metadata(path)
        tile["frame_index"] = index
        if tile["channels"] != 1 or tile["frames"] != 1:
            raise ValueError(
                f"{path.name}: one channel and time point required"
            )
        tiles.append(tile)
    first = tiles[0]
    for tile in tiles:
        for axis in ("x", "y"):
            key = f"pixel_size_{axis}_um"
            size = tile[key]
            if not math.isfinite(size) or size <= 0:
                raise ValueError(f"{tile['filename']}: invalid {key}")
            if not math.isclose(size, first[key], rel_tol=1e-6, abs_tol=1e-9):
                raise ValueError(
                    f"Inconsistent {key} across the nine original tiles"
                )
        for key in ("width_px", "height_px", "slices"):
            if tile[key] <= 0 or tile[key] != first[key]:
                raise ValueError(f"Inconsistent tile dimensions: {key}")
    x_size = first["pixel_size_x_um"]
    y_size = first["pixel_size_y_um"]
    return {
        "pixel_size_x_um": x_size,
        "pixel_size_y_um": y_size,
        "pixel_area_um2": x_size * y_size,
        "tiles": tiles,
    }
