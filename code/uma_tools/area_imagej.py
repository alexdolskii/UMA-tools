"""
ImageJ operations for native-resolution FN projection and area masks.
"""

from __future__ import annotations

import importlib.metadata
import json
import math
import struct

from .files import sha256_file
from .imagej import (
    FIJI_ENDPOINT,
    ImageJInitializationError,
    initialize_imagej,
)

FLOAT32_MAX = 3.4028234663852886e38
DEFAULT_THRESHOLD_LOWER = 2000.0


class ValidationError(Exception):
    """
    An input or selection problem that must not be silently bypassed.
    """


def float32_limit(value, label):
    """
    Report the effective precision used by ImageJ's 32-bit thresholding.
    """
    try:
        number = float(value)
        if not math.isfinite(number) or number < 0 or number > FLOAT32_MAX:
            raise ValueError("outside the finite non-negative float32 range")
        effective = struct.unpack("!f", struct.pack("!f", number))[0]
    except (ValueError, TypeError, OverflowError, struct.error) as error:
        raise ValidationError(
            f"{label} must be finite, non-negative, "
            "and representable in 32-bit float."
        ) from error
    return number, effective


def threshold_settings(lower, upper):
    if lower is None:
        if upper is not None:
            raise ValidationError(
                "An upper threshold without a lower threshold is not valid."
            )
        return None
    requested_lower, effective_lower = float32_limit(lower, "Lower threshold")
    unbounded = upper is None or str(upper).strip().lower() in (
        "inf",
        "+inf",
        "infinity",
        "none",
    )
    if unbounded:
        requested_upper, effective_upper = None, FLOAT32_MAX
    else:
        requested_upper, effective_upper = float32_limit(
            upper, "Upper threshold"
        )
        if requested_upper < requested_lower:
            raise ValidationError(
                "Upper threshold must be greater than or equal "
                "to the lower threshold."
            )
    return {
        "requested_lower": requested_lower,
        "requested_upper": requested_upper,
        "lower": effective_lower,
        "upper": effective_upper,
        "upper_unbounded": unbounded,
    }


class ImageJEngine:
    """Read, project, and threshold original fields as float pixels."""

    def __init__(self, log):
        log.event(
            "INFO",
            "ImageJ initialization",
            f"Starting imagej.init('{FIJI_ENDPOINT}', mode='headless').",
        )
        self.ij = initialize_imagej()
        try:
            import numpy as np
            import scyjava as sj

            self.np = np
            self.IJ = sj.jimport("ij.IJ")
            self.ImagePlus = sj.jimport("ij.ImagePlus")
            self.ImageStack = sj.jimport("ij.ImageStack")
            self.ZProjector = sj.jimport("ij.plugin.ZProjector")
            self.FileSaver = sj.jimport("ij.io.FileSaver")
            self.ImageProcessor = sj.jimport("ij.process.ImageProcessor")
            self.Measurements = sj.jimport("ij.measure.Measurements")
            self.sj = sj
            self.BF = self.ImporterOptions = self.ImageReader = None
            system = sj.jimport("java.lang.System")
            self.versions = {
                "ImageJ1": str(self.IJ.getVersion()),
                "Java": str(system.getProperty("java.version")),
                "pyimagej": importlib.metadata.version("pyimagej"),
                "scyjava": importlib.metadata.version("scyjava"),
                "numpy": np.__version__,
                "imagej_endpoint": FIJI_ENDPOINT,
                "imagej_mode": "headless",
            }
        except Exception:
            self.close()
            raise
        log.event("INFO", "ImageJ versions", json.dumps(self.versions))

    def load_bioformats(self):
        if self.BF is None:
            try:
                self.BF = self.sj.jimport("loci.plugins.BF")
                self.ImporterOptions = self.sj.jimport(
                    "loci.plugins.in.ImporterOptions"
                )
                self.ImageReader = self.sj.jimport("loci.formats.ImageReader")
                tools = self.sj.jimport("loci.formats.FormatTools")
                self.versions["Bio-Formats"] = str(tools.VERSION)
            except Exception as error:
                raise ImageJInitializationError(
                    "Fiji's Bio-Formats classes are unavailable "
                    "for ND2 import: " + str(error)
                ) from error

    @staticmethod
    def release(imp):
        if imp is not None:
            imp.changes = False
            imp.close()

    def open_original(self, path, channel):
        if path.suffix.lower() == ".nd2":
            self.load_bioformats()
            reader = self.ImageReader()
            try:
                reader.setGroupFiles(False)
                reader.setId(str(path))
                if int(reader.getSeriesCount()) != 1:
                    raise ValidationError(
                        f"ND2 contains {reader.getSeriesCount()} series: "
                        f"{path.name}. Export individual fields first."
                    )
                reader.setSeries(0)
                if int(reader.getSizeT()) != 1:
                    raise ValidationError(
                        "Multiple time points are not supported "
                        f"for one image row: {path.name}"
                    )
                if bool(reader.isRGB()):
                    raise ValidationError(
                        "RGB-packed source data are not supported: "
                        f"{path.name}"
                    )
                if channel > int(reader.getSizeC()):
                    raise ValidationError(
                        f"FN channel {channel} exceeds "
                        f"{reader.getSizeC()} channels in {path.name}"
                    )
                expected = [
                    int(reader.getSizeX()),
                    int(reader.getSizeY()),
                    int(reader.getSizeC()),
                    int(reader.getSizeZ()),
                    int(reader.getSizeT()),
                ]
            finally:
                reader.close()
            options = self.ImporterOptions()
            options.setId(str(path))
            options.setAutoscale(False)
            options.setOpenAllSeries(False)
            options.setQuiet(True)
            options.setWindowless(True)
            options.setColorMode(self.ImporterOptions.COLOR_MODE_GRAYSCALE)
            # ImporterOptions inherits ImageJ preferences. Explicitly
            # disable
            # saved cropping, splitting and range settings for
            # reproducibility.
            options.setGroupFiles(False)
            options.setCrop(False)
            options.setSpecifyRanges(False)
            options.setSplitChannels(False)
            options.setSplitFocalPlanes(False)
            options.setSplitTimepoints(False)
            options.setVirtual(False)
            options.setSwapDimensions(False)
            options.setConcatenate(False)
            options.setShowMetadata(False)
            options.setShowOMEXML(False)
            options.setShowROIs(False)
            options.setStackFormat(self.ImporterOptions.VIEW_HYPERSTACK)
            options.setStackOrder(self.ImporterOptions.ORDER_XYCZT)
            options.clearSeries()
            options.setSeriesOn(0, True)
            images = self.BF.openImagePlus(options)
            if images is None or len(images) != 1:
                if images is not None:
                    for image in images:
                        self.release(image)
                raise ValidationError(
                    f"Expected exactly one imported field from {path.name}"
                )
            imp = images[0]
            if [int(v) for v in imp.getDimensions()] != expected:
                self.release(imp)
                raise ValidationError(
                    "Bio-Formats metadata and imported dimensions "
                    f"disagree: {path.name}"
                )
            return imp, "Bio-Formats ND2 (autoscale disabled)"
        imp = self.IJ.openImage(str(path))
        if imp is None:
            raise ValidationError(
                f"ImageJ could not open original TIFF: {path}"
            )
        return imp, "ImageJ TIFF"

    def project(self, imp, channel, method):
        width, height, channels, slices, frames = [
            int(v) for v in imp.getDimensions()
        ]
        if width <= 0 or height <= 0 or slices < 1 or frames != 1:
            raise ValidationError(
                "Expected one field with one time point "
                "and at least one Z slice; "
                f"got {list(imp.getDimensions())}"
            )
        if not 1 <= channel <= channels or int(imp.getBitDepth()) not in (
            8,
            16,
            32,
        ):
            raise ValidationError(
                f"Invalid channel {channel} or unsupported source "
                f"bit depth {imp.getBitDepth()}."
            )
        if int(imp.getStackSize()) != channels * slices:
            raise ValidationError(
                "Original stack dimensions do not match "
                "C * Z for one time point."
            )
        calibration = imp.getCalibration().copy()
        if bool(calibration.calibrated()):
            raise ValidationError(
                "An intensity-calibrated source requires an explicit "
                "intensity-unit policy. Raw uncalibrated intensities "
                "are required."
            )
        px, py = float(calibration.pixelWidth), float(calibration.pixelHeight)
        if not all(math.isfinite(v) and v > 0 for v in (px, py)):
            raise ValidationError(
                "Original image has invalid spatial calibration."
            )
        stack = self.ImageStack(width, height)
        for z in range(1, slices + 1):
            processor = (
                imp.getStack()
                .getProcessor(imp.getStackIndex(channel, z, 1))
                .convertToFloatProcessor()
            )
            values = self.np.asarray(
                processor.getPixels(), dtype=self.np.float32
            )
            if not self.np.isfinite(values).all():
                raise ValidationError(
                    "Non-finite input intensities "
                    f"in channel {channel}, slice {z}."
                )
            stack.addSlice(processor)
        selected = self.ImagePlus("FN selected channel", stack)
        selected.setDimensions(1, slices, 1)
        selected.setCalibration(calibration.copy())
        try:
            if slices == 1:
                projection = self.ImagePlus(
                    "FN projection", stack.getProcessor(1).duplicate()
                )
            else:
                projector = self.ZProjector(selected)
                modes = {
                    "sum": self.ZProjector.SUM_METHOD,
                    "mean": self.ZProjector.AVG_METHOD,
                    "max": self.ZProjector.MAX_METHOD,
                }
                projector.setMethod(modes[method])
                projector.setStartSlice(1)
                projector.setStopSlice(slices)
                projector.doProjection()
                projection = projector.getProjection()
            if projection is None or int(projection.getBitDepth()) != 32:
                self.release(projection)
                raise RuntimeError(
                    "ImageJ did not produce a 32-bit projection."
                )
            projection.setCalibration(calibration.copy())
            projection.deleteRoi()
            processor = projection.getProcessor()
            processor.resetThreshold()
            values = self.np.asarray(
                processor.getPixels(), dtype=self.np.float32
            )
            if not self.np.isfinite(values).all():
                self.release(projection)
                raise ValidationError(
                    "Non-finite intensities or float32 overflow "
                    "in the projection."
                )
            unit = str(calibration.getUnit()) or "pixel"
            metadata = {
                "Width_Pixels": width,
                "Height_Pixels": height,
                "Bit_Depth": 32,
                "Source_Bit_Depth": int(imp.getBitDepth()),
                "Source_Channels": channels,
                "Channel_Index": channel,
                "Number_of_Z_Stacks": slices,
                "Source_Timepoints": 1,
                "Projection_Method": method.upper(),
                "Projection_Min": float(values.min()),
                "Projection_Max": float(values.max()),
                "Total_Pixels": width * height,
                "Image_Area": width * height * px * py,
                "Area_Unit": unit + "^2",
                "Pixel_Width": px,
                "Pixel_Height": py,
                "Pixel_Unit": unit,
            }
            return projection, metadata
        finally:
            self.release(selected)

    def save_projection(self, projection, path):
        """Verify the saved TIFF values, dimensions, and calibration."""
        if (
            not self.FileSaver(projection).saveAsTiff(str(path))
            or not path.is_file()
        ):
            raise RuntimeError(f"Failed to save 32-bit projection: {path}")
        reopened = self.IJ.openImage(str(path))
        try:
            if (
                reopened is None
                or int(reopened.getBitDepth()) != 32
                or [int(v) for v in reopened.getDimensions()]
                != [int(v) for v in projection.getDimensions()]
            ):
                raise RuntimeError(
                    "The stored projection is not the expected "
                    f"single-plane 32-bit TIFF: {path}"
                )
            expected = self.np.asarray(
                projection.getProcessor().getPixels(), dtype=self.np.float32
            )
            actual = self.np.asarray(
                reopened.getProcessor().getPixels(), dtype=self.np.float32
            )
            if not self.np.array_equal(expected, actual):
                raise RuntimeError(
                    "TIFF export changed numeric projection "
                    f"intensities: {path}"
                )
            a, b = projection.getCalibration(), reopened.getCalibration()
            if (
                not math.isclose(
                    float(a.pixelWidth), float(b.pixelWidth), rel_tol=1e-6
                )
                or not math.isclose(
                    float(a.pixelHeight), float(b.pixelHeight), rel_tol=1e-6
                )
                or str(a.getUnit()) != str(b.getUnit())
            ):
                raise RuntimeError(
                    f"TIFF export changed spatial calibration: {path}"
                )
        finally:
            self.release(reopened)
        return sha256_file(path)

    def measure(self, projection, mask_path, limits):
        """Threshold float pixels without 8-bit intensity conversion."""
        processor = projection.getProcessor()
        processor.setThreshold(
            limits["lower"], limits["upper"], self.ImageProcessor.NO_LUT_UPDATE
        )
        values = self.np.asarray(processor.getPixels(), dtype=self.np.float32)
        expected = (values >= limits["lower"]) & (values <= limits["upper"])
        mask_processor = processor.createMask()
        if mask_processor is None:
            raise RuntimeError(
                "ImageJ could not create a mask "
                "from the explicit float threshold."
            )
        mask = self.ImagePlus(mask_path.name, mask_processor)
        mask.setCalibration(projection.getCalibration().copy())
        try:
            mask_values = self.np.asarray(
                mask_processor.getPixels(), dtype=self.np.int8
            ).view(self.np.uint8)
            if not self.np.array_equal(
                mask_values, expected.astype(self.np.uint8) * 255
            ):
                raise RuntimeError(
                    "ImageJ float-threshold mask disagrees "
                    "with the raw projection values."
                )
            mask_processor.setThreshold(
                255.0, 255.0, self.ImageProcessor.NO_LUT_UPDATE
            )
            flags = (
                int(self.Measurements.AREA)
                | int(self.Measurements.LIMIT)
                | int(self.Measurements.AREA_FRACTION)
            )
            stats = mask.getStatistics(flags)
            positive = int(stats.pixelCount)
            total = int(values.size)
            percentage = float(stats.areaFraction)
            area = float(stats.area)
            cal = projection.getCalibration()
            if (
                positive != int(expected.sum())
                or not math.isclose(
                    percentage,
                    positive / total * 100,
                    abs_tol=1e-10,
                    rel_tol=1e-10,
                )
                or not math.isclose(
                    area,
                    positive * float(cal.pixelWidth) * float(cal.pixelHeight),
                    abs_tol=1e-9,
                    rel_tol=1e-9,
                )
            ):
                raise RuntimeError(
                    "ImageJ area measurements failed reconciliation."
                )
            mask_processor.resetThreshold()
            if not self.FileSaver(mask).saveAsTiff(str(mask_path)):
                raise RuntimeError(f"Failed to save binary mask: {mask_path}")
            reopened = self.IJ.openImage(str(mask_path))
            try:
                if (
                    reopened is None
                    or int(reopened.getBitDepth()) != 8
                    or not self.np.array_equal(
                        self.np.asarray(
                            reopened.getProcessor().getPixels(),
                            dtype=self.np.int8,
                        ).view(self.np.uint8),
                        mask_values,
                    )
                ):
                    raise RuntimeError(
                        "Saved TIFF mask pixels differ "
                        "from measured mask pixels."
                    )
            finally:
                self.release(reopened)
            return {
                "FN_Positive_Pixels": positive,
                "FN_Area_Percent": percentage,
                "FN_Area": area,
                "Threshold_Lower": limits["lower"],
                "Threshold_Upper": None
                if limits["upper_unbounded"]
                else limits["upper"],
                "Effective_Threshold_Upper": limits["upper"],
                "Threshold_Units": "Raw projection intensity",
                "Requested_Threshold_Lower": limits["requested_lower"],
                "Requested_Threshold_Upper": limits["requested_upper"],
                "Mask_Bit_Depth": 8,
                "Mask_File_Name": mask_path.name,
                "Mask_Path": str(mask_path),
                "Mask_SHA256": sha256_file(mask_path),
            }
        finally:
            processor.resetThreshold()
            self.release(mask)

    def close(self):
        if self.ij is not None:
            try:
                self.ij.dispose()
            finally:
                self.ij = None
