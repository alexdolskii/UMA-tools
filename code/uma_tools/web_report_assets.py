"""Bounded-memory PNG previews of saved two-dimensional images only."""

from __future__ import annotations

import numpy as np
from PIL import Image

from .files import sha256_file
from .web_report_data import Resolver, uid


class Assets:
    def __init__(self, source, output, logger):
        self.resolver = Resolver(source)
        self.output = output
        self.logger = logger
        self.items = {}
        self.ranges = {}
        self.failures = []

    def add(self, value, base, label, cohort, mode="linear", expected=None):
        path = self.resolver.resolve(value, base)
        source = str(value)
        key = uid(f"{path or source}|{base}|{cohort}|{mode}|{expected}")
        if key not in self.items:
            item = {
                "id": key,
                "source": str(path or source),
                "path": path,
                "label": label,
                "cohort": cohort,
                "mode": mode,
                "expected_sha256": expected,
                "status": "PENDING",
            }
            self.items[key] = item
            if path is None:
                self.fail(item, "Recorded file missing, unsafe or ambiguous")
        return {"asset": key, "label": label, "source": str(path or source)}

    def fail(self, item, message):
        item.update(status="UNAVAILABLE", error=str(message))
        self.failures.append(f"{item['source']}: {message}")
        self.logger.warning(
            "PREVIEW_UNAVAILABLE | %s | %s", item["source"], message
        )

    @staticmethod
    def array(item):
        path = item["path"]
        if path.suffix.lower() not in (".tif", ".tiff", ".png"):
            raise ValueError("Only saved TIFF/PNG images can be previewed")
        with Image.open(path) as image:
            if getattr(image, "n_frames", 1) != 1:
                raise ValueError("A stack is not a saved 2D projection")
            if image.mode == "P":
                image = image.convert("RGB")
            array = np.array(image)
        if array.ndim != 2 and not (
            array.ndim == 3 and array.shape[-1] in (3, 4)
        ):
            raise ValueError(f"Expected 2D grayscale/RGB; got {array.shape}")
        return array

    def build(self, progress):
        folder = self.output / "assets" / "previews"
        folder.mkdir(parents=True)
        for index, item in enumerate(self.items.values(), 1):
            progress(index, len(self.items), "Inspecting saved previews")
            if item["status"] != "PENDING":
                continue
            try:
                checksum = sha256_file(item["path"])
                if (
                    item["expected_sha256"]
                    and checksum != (item["expected_sha256"])
                ):
                    raise ValueError("Recorded image SHA256 does not match")
                array = self.array(item)
                key = item["cohort"] + "|" + str(array.dtype)
                item.update(
                    sha256=checksum,
                    width=array.shape[1],
                    height=array.shape[0],
                    dtype=str(array.dtype),
                    range_key=key,
                    status="READY",
                )
                if array.ndim == 2 and item["mode"] == "linear":
                    finite = array[np.isfinite(array)]
                    if not finite.size:
                        raise ValueError("No finite pixels in saved image")
                    low, high = float(finite.min()), float(finite.max())
                    previous = self.ranges.get(key, [low, high])
                    self.ranges[key] = [
                        min(low, previous[0]),
                        max(high, previous[1]),
                    ]
            except (
                OSError,
                ValueError,
                Image.DecompressionBombError,
            ) as error:
                self.fail(item, error)
                self.logger.exception("PREVIEW_READ_FAILED")
        for index, item in enumerate(self.items.values(), 1):
            progress(index, len(self.items), "Saving PNG previews")
            if item["status"] != "READY":
                continue
            try:
                if sha256_file(item["path"]) != item["sha256"]:
                    raise ValueError("Image changed during preview generation")
                array = self.array(item)
                if array.ndim == 3:
                    if array.dtype != np.uint8:
                        raise ValueError("RGB preview must be uint8")
                elif item["mode"] == "binary":
                    array = np.where(array != 0, 255, 0).astype(np.uint8)
                elif item["mode"] == "linear":
                    low, high = self.ranges[item["range_key"]]
                    item["display_range"] = [low, high]
                    values = np.nan_to_num(
                        array.astype(float), nan=low, posinf=high, neginf=low
                    )
                    array = (
                        np.zeros(array.shape, np.uint8)
                        if high == low
                        else (
                            np.clip(
                                (values - low) * 255 / (high - low), 0, 255
                            ).astype(np.uint8)
                        )
                    )
                image = Image.fromarray(array)
                image.thumbnail((1024, 1024), Image.Resampling.NEAREST)
                target = folder / (item["id"] + ".png")
                image.save(target)
                if sha256_file(item["path"]) != item["sha256"]:
                    target.unlink()
                    raise ValueError("Image changed while saving preview")
                item.update(
                    preview=target.relative_to(self.output).as_posix(),
                    preview_width=image.width,
                    preview_height=image.height,
                    status="AVAILABLE",
                    preview_sha256=sha256_file(target),
                )
                self.logger.info(
                    "PREVIEW_SAVED | %s | range=%s",
                    item["source"],
                    item.get("display_range"),
                )
            except (
                OSError,
                ValueError,
                Image.DecompressionBombError,
            ) as error:
                self.fail(item, error)
                self.logger.exception("PREVIEW_SAVE_FAILED")
        return {
            key: {k: v for k, v in item.items() if k != "path"}
            for key, item in self.items.items()
        }
