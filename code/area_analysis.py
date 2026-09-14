#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Measure FN area from original ND2/TIFF images using native SUM32 projections.

Run in the UMA environment:
    area_analysis -i input_paths.json
    area_analysis -i input_paths.json -t 2000
    area_analysis -i input_paths.json -t 2000 50000

The threshold is an inclusive LOWER [UPPER] interval in raw projection intensity
units. The default lower bound is 2000; an omitted upper bound (or inf) uses the
largest finite float32 value. The 1-based fibronectin channel is requested once
unless --channel is supplied. --projection defaults to sum; --projections-only
retains the optional inspection mode without masks or area measurements.

Every source folder in the folder_paths JSON is scanned directly for visible
.nd2, .tif and .tiff files. Alignment/thickness outputs and _Seq identifiers are
not required. Subfolders, including previous results, are not scanned.
Image_ID is the complete original filename, including its extension.

Each source folder receives a unique Area_assay_results_<timestamp>_<id> folder
with SUM32 projections, 0/255 masks, area tables, parameters, status and logs.
Startup failures use a new results folder in an available source folder, or in
the working directory if no source folder is available. No logs are written to
the installed package. An image error stops its folder and retains partial
results; other source folders continue. Command errors return a nonzero status.

Original float intensities, native XY resolution, spatial calibration, SUM
projection and full-frame area measurements are preserved. No resizing,
intensity scaling, denoising or background subtraction is applied. Multi-series
ND2, RGB-packed data and multiple time points remain unsupported. SUM includes
background and depends on Z count; the requested/effective thresholds and Z
count are recorded for every image. Fiji uses sc.fiji:fiji:2.14.0 in headless
mode, matching the other UMA commands.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import struct
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_VERSION = "2.1.0"
FIJI_ENDPOINT = "sc.fiji:fiji:2.14.0"
DEFAULT_THRESHOLD_LOWER = 2000.0
ORIGINAL_EXTENSIONS = (".nd2", ".tif", ".tiff")
EVENT_COLUMNS = ["Timestamp_UTC", "Level", "Stage", "Message"]
FLOAT32_MAX = 3.4028234663852886e38
MANIFEST_COLUMNS = ["File_Name", "Image_ID", "Selected", "Reason", "Path", "Bytes", "SHA256"]
PROJECTION_COLUMNS = [
    "File_Name", "Image_ID", "Width_Pixels", "Height_Pixels", "Bit_Depth",
    "Source_Bit_Depth", "Source_Channels", "Channel_Index", "Number_of_Z_Stacks",
    "Source_Timepoints", "Projection_Method", "Projection_Min", "Projection_Max",
    "Total_Pixels", "Image_Area", "Area_Unit", "Pixel_Width", "Pixel_Height", "Pixel_Unit",
    "Source_Original_Path", "Source_Projection_Path", "Projection_SHA256",
    "Source_Folder", "Source_Reader", "Source_SHA256", "Program_Version", "Run_ID"
]
SUMMARY_COLUMNS = PROJECTION_COLUMNS + [
    "FN_Positive_Pixels", "FN_Area_Percent", "FN_Area", "Threshold_Lower", "Threshold_Upper",
    "Effective_Threshold_Upper", "Threshold_Units", "Requested_Threshold_Lower", "Requested_Threshold_Upper",
    "Mask_Bit_Depth", "Mask_File_Name", "Mask_Path", "Mask_SHA256"
]


class ValidationError(Exception):
    """An input or selection problem that must not be silently bypassed."""


class ImageJInitializationError(Exception):
    """Exception raised for unsuccessful initialization of ImageJ."""


def package_version():
    try:
        return importlib.metadata.version("uma-tools")
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def initialize_imagej():
    """Initialize the same headless Fiji endpoint as alignment and thickness."""
    print("Initializing ImageJ...", flush=True)
    try:
        import imagej

        ij = imagej.init(FIJI_ENDPOINT, mode="headless")
    except Exception as error:
        raise ImageJInitializationError(f"Failed to initialize ImageJ: {error}") from error
    print("ImageJ initialization completed.", flush=True)
    return ij


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def new_output_folder(parent):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    base = f"{stamp}_{os.getpid()}"
    for counter in range(10000):
        run_id = base + (f"_{counter:03d}" if counter else "")
        path = parent / f"Area_assay_results_{run_id}"
        try:
            path.mkdir()
            return run_id, path
        except FileExistsError:
            continue
    raise RuntimeError("Could not allocate a unique results folder.")


def save_json(path, value):
    temporary = path.with_name(path.name + ".pending")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def save_csv(path, columns, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class RunLog:
    def __init__(self, directory, append=False):
        mode = "a" if append else "w"
        csv_path = directory / "run_log.csv"
        needs_header = not append or not csv_path.exists() or csv_path.stat().st_size == 0
        self.text_stream = (directory / "run.log").open(mode, encoding="utf-8")
        self.csv_stream = csv_path.open(mode, encoding="utf-8-sig", newline="")
        self.writer = csv.DictWriter(self.csv_stream, fieldnames=EVENT_COLUMNS)
        if needs_header:
            self.writer.writeheader()
        self.csv_stream.flush()

    def event(self, level, stage, message):
        record = dict(zip(EVENT_COLUMNS, [utc_now(), level, stage, str(message)]))
        line = f"[{record['Timestamp_UTC']}] [{level}] [{stage}] {message}"
        self.text_stream.write(line + "\n")
        self.text_stream.flush()
        self.writer.writerow(record)
        self.csv_stream.flush()
        print(line, flush=True)

    def close(self):
        self.text_stream.close()
        self.csv_stream.close()


def resolve_path(value, relative_to):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else relative_to / path).resolve()


def read_source_folders(args):
    """Read explicit source paths without depending on the installed package."""
    if args.folder is not None:
        folders = [resolve_path(args.folder, Path.cwd())]
        input_json = None
    else:
        if Path(args.input).name.startswith("._"):
            raise ValidationError(f"macOS metadata files cannot be used as input: {args.input}")
        input_json = resolve_path(args.input, Path.cwd())
        if input_json.name.startswith("._"):
            raise ValidationError(f"macOS metadata files cannot be used as input: {input_json}")
        try:
            with input_json.open(encoding="utf-8-sig") as stream:
                value = json.load(stream)
        except FileNotFoundError as error:
            raise ValidationError(f"Input JSON was not found: {input_json}") from error
        values = value.get("folder_paths") if isinstance(value, dict) else None
        if not isinstance(values, list) or not values:
            raise ValidationError("The JSON must contain a nonempty folder_paths list.")
        if any(not isinstance(item, str) or not item.strip() for item in values):
            raise ValidationError("Every folder_paths entry must be a nonempty path string.")
        folders = [resolve_path(item, input_json.parent) for item in values]
    if len(set(folders)) != len(folders):
        raise ValidationError("The input repeats the same source folder. No folder was silently removed.")
    return folders, input_json


def original_inventory(source_folder):
    """Select every visible original image directly in the supplied folder."""
    files = sorted(path for path in source_folder.iterdir()
                   if not path.name.startswith(".") and path.is_file()
                   and path.suffix.lower() in ORIGINAL_EXTENSIONS)
    if not files:
        raise ValidationError(f"No visible ND2/TIF/TIFF images found in: {source_folder}")
    return [{"File_Name": path.name, "Image_ID": path.name, "Selected": True,
             "Reason": "DIRECT_SOURCE_IMAGE", "Path": str(path),
             "Bytes": path.stat().st_size, "SHA256": ""} for path in files]


def save_startup_error(error, folders, args):
    """Create a results directory for failures before image processing starts."""
    parents = [folder for folder in folders if folder.is_dir()]
    parents.append(Path.cwd())
    for parent in dict.fromkeys(parents):
        try:
            run_id, output = new_output_folder(parent)
        except OSError:
            continue
        log = None
        try:
            log = RunLog(output)
            log.event("ERROR", "Startup", str(error) or "Run interrupted.")
            (output / "traceback.txt").write_text(traceback.format_exc(), encoding="utf-8")
            save_csv(output / "errors.csv", ["File_Name", "Stage", "Issue"],
                     [{"File_Name": "", "Stage": "Startup", "Issue": str(error) or "Run interrupted."}])
            state = "CANCELLED" if isinstance(error, KeyboardInterrupt) else "ERROR"
            save_json(output / "run_status.json",
                      {"run_id": run_id, "script_version": SCRIPT_VERSION,
                       "package_version": package_version(), "status": state,
                       "stage": "Startup", "started_utc": utc_now(), "ended_utc": utc_now(),
                       "processed_images": 0, "error": str(error) or "Run interrupted.",
                       "run_directory": str(output)})
            save_json(output / "run_parameters.json",
                      {"run_id": run_id, "script_version": SCRIPT_VERSION,
                       "package_version": package_version(), "arguments": vars(args),
                       "python": platform.python_version(), "platform": platform.platform(),
                       "python_executable": sys.executable, "imagej_endpoint": FIJI_ENDPOINT})
            print(f"Startup diagnostics: {output}", file=sys.stderr, flush=True)
            return output
        except OSError as log_error:
            print(f"Could not write startup diagnostics in {output}: {log_error}",
                  file=sys.stderr, flush=True)
        finally:
            if log is not None:
                log.close()
    print("Could not create a writable results directory for startup diagnostics.",
          file=sys.stderr, flush=True)
    return None


def float32_limit(value, label):
    """Report the effective precision used by ImageJ's 32-bit thresholding."""
    try:
        number = float(value)
        if not math.isfinite(number) or number < 0 or number > FLOAT32_MAX:
            raise ValueError("outside the finite non-negative float32 range")
        effective = struct.unpack('!f', struct.pack('!f', number))[0]
    except (ValueError, TypeError, OverflowError, struct.error) as error:
        raise ValidationError(f"{label} must be finite, non-negative, and representable in 32-bit float.") from error
    return number, effective


def threshold_settings(lower, upper):
    if lower is None:
        if upper is not None:
            raise ValidationError("An upper threshold without a lower threshold is not valid.")
        return None
    requested_lower, effective_lower = float32_limit(lower, "Lower threshold")
    unbounded = upper is None or str(upper).strip().lower() in ("inf", "+inf", "infinity", "none")
    if unbounded:
        requested_upper, effective_upper = None, FLOAT32_MAX
    else:
        requested_upper, effective_upper = float32_limit(upper, "Upper threshold")
        if requested_upper < requested_lower:
            raise ValidationError("Upper threshold must be greater than or equal to the lower threshold.")
    return {"requested_lower": requested_lower, "requested_upper": requested_upper,
            "lower": effective_lower, "upper": effective_upper, "upper_unbounded": unbounded}


class ImageJEngine:
    """Read original fields, project the chosen channel, and threshold float pixels."""

    def __init__(self, log):
        log.event("INFO", "ImageJ initialization", f"Starting imagej.init('{FIJI_ENDPOINT}', mode='headless').")
        self.ij = initialize_imagej()
        try:
            import numpy as np
            import scyjava as sj
            self.np = np
            self.IJ = sj.jimport('ij.IJ')
            self.ImagePlus = sj.jimport('ij.ImagePlus')
            self.ImageStack = sj.jimport('ij.ImageStack')
            self.ZProjector = sj.jimport('ij.plugin.ZProjector')
            self.FileSaver = sj.jimport('ij.io.FileSaver')
            self.ImageProcessor = sj.jimport('ij.process.ImageProcessor')
            self.Measurements = sj.jimport('ij.measure.Measurements')
            self.sj = sj
            self.BF = self.ImporterOptions = self.ImageReader = None
            system = sj.jimport('java.lang.System')
            self.versions = {"ImageJ1": str(self.IJ.getVersion()), "Java": str(system.getProperty('java.version')),
                             "pyimagej": importlib.metadata.version('pyimagej'),
                             "scyjava": importlib.metadata.version('scyjava'), "numpy": np.__version__,
                             "imagej_endpoint": FIJI_ENDPOINT, "imagej_mode": "headless"}
        except Exception:
            self.close()
            raise
        log.event("INFO", "ImageJ versions", json.dumps(self.versions))

    def load_bioformats(self):
        if self.BF is None:
            try:
                self.BF = self.sj.jimport('loci.plugins.BF')
                self.ImporterOptions = self.sj.jimport('loci.plugins.in.ImporterOptions')
                self.ImageReader = self.sj.jimport('loci.formats.ImageReader')
                tools = self.sj.jimport('loci.formats.FormatTools')
                self.versions['Bio-Formats'] = str(tools.VERSION)
            except Exception as error:
                raise ImageJInitializationError("Fiji's Bio-Formats classes are unavailable for ND2 import: " + str(error)) from error

    @staticmethod
    def release(imp):
        if imp is not None:
            imp.changes = False
            imp.close()

    def open_original(self, path, channel):
        if path.suffix.lower() == '.nd2':
            self.load_bioformats()
            reader = self.ImageReader()
            try:
                reader.setGroupFiles(False)
                reader.setId(str(path))
                if int(reader.getSeriesCount()) != 1:
                    raise ValidationError(f"ND2 contains {reader.getSeriesCount()} series: {path.name}. Export individual fields first.")
                reader.setSeries(0)
                if int(reader.getSizeT()) != 1:
                    raise ValidationError(f"Multiple time points are not supported for one image row: {path.name}")
                if bool(reader.isRGB()):
                    raise ValidationError(f"RGB-packed source data are not supported: {path.name}")
                if channel > int(reader.getSizeC()):
                    raise ValidationError(f"FN channel {channel} exceeds {reader.getSizeC()} channels in {path.name}")
                expected = [int(reader.getSizeX()), int(reader.getSizeY()), int(reader.getSizeC()),
                            int(reader.getSizeZ()), int(reader.getSizeT())]
            finally:
                reader.close()
            options = self.ImporterOptions()
            options.setId(str(path))
            options.setAutoscale(False)
            options.setOpenAllSeries(False)
            options.setQuiet(True)
            options.setWindowless(True)
            options.setColorMode(self.ImporterOptions.COLOR_MODE_GRAYSCALE)
            # ImporterOptions inherits ImageJ preferences. Explicitly disable
            # saved cropping, splitting and range settings for reproducibility.
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
                raise ValidationError(f"Expected exactly one imported field from {path.name}")
            imp = images[0]
            if [int(v) for v in imp.getDimensions()] != expected:
                self.release(imp)
                raise ValidationError(f"Bio-Formats metadata and imported dimensions disagree: {path.name}")
            return imp, 'Bio-Formats ND2 (autoscale disabled)'
        imp = self.IJ.openImage(str(path))
        if imp is None:
            raise ValidationError(f"ImageJ could not open original TIFF: {path}")
        return imp, 'ImageJ TIFF'

    def project(self, imp, channel, method):
        width, height, channels, slices, frames = [int(v) for v in imp.getDimensions()]
        if width <= 0 or height <= 0 or slices < 1 or frames != 1:
            raise ValidationError(f"Expected one field with one time point and at least one Z slice; got {list(imp.getDimensions())}")
        if not 1 <= channel <= channels or int(imp.getBitDepth()) not in (8, 16, 32):
            raise ValidationError(f"Invalid channel {channel} or unsupported source bit depth {imp.getBitDepth()}.")
        if int(imp.getStackSize()) != channels * slices:
            raise ValidationError("Original stack dimensions do not match C * Z for one time point.")
        calibration = imp.getCalibration().copy()
        if bool(calibration.calibrated()):
            raise ValidationError("An intensity-calibrated source requires an explicit intensity-unit policy. Raw uncalibrated intensities are required.")
        px, py = float(calibration.pixelWidth), float(calibration.pixelHeight)
        if not all(math.isfinite(v) and v > 0 for v in (px, py)):
            raise ValidationError("Original image has invalid spatial calibration.")
        stack = self.ImageStack(width, height)
        for z in range(1, slices + 1):
            processor = imp.getStack().getProcessor(imp.getStackIndex(channel, z, 1)).convertToFloatProcessor()
            values = self.np.asarray(processor.getPixels(), dtype=self.np.float32)
            if not self.np.isfinite(values).all():
                raise ValidationError(f"Non-finite input intensities in channel {channel}, slice {z}.")
            stack.addSlice(processor)
        selected = self.ImagePlus('FN selected channel', stack)
        selected.setDimensions(1, slices, 1)
        selected.setCalibration(calibration.copy())
        try:
            if slices == 1:
                projection = self.ImagePlus('FN projection', stack.getProcessor(1).duplicate())
            else:
                projector = self.ZProjector(selected)
                modes = {'sum': self.ZProjector.SUM_METHOD, 'mean': self.ZProjector.AVG_METHOD,
                         'max': self.ZProjector.MAX_METHOD}
                projector.setMethod(modes[method])
                projector.setStartSlice(1)
                projector.setStopSlice(slices)
                projector.doProjection()
                projection = projector.getProjection()
            if projection is None or int(projection.getBitDepth()) != 32:
                self.release(projection)
                raise RuntimeError("ImageJ did not produce a 32-bit projection.")
            projection.setCalibration(calibration.copy())
            projection.deleteRoi()
            processor = projection.getProcessor()
            processor.resetThreshold()
            values = self.np.asarray(processor.getPixels(), dtype=self.np.float32)
            if not self.np.isfinite(values).all():
                self.release(projection)
                raise ValidationError("Non-finite intensities or float32 overflow in the projection.")
            unit = str(calibration.getUnit()) or 'pixel'
            metadata = {"Width_Pixels": width, "Height_Pixels": height, "Bit_Depth": 32,
                        "Source_Bit_Depth": int(imp.getBitDepth()), "Source_Channels": channels,
                        "Channel_Index": channel, "Number_of_Z_Stacks": slices, "Source_Timepoints": 1,
                        "Projection_Method": method.upper(), "Projection_Min": float(values.min()),
                        "Projection_Max": float(values.max()), "Total_Pixels": width * height,
                        "Image_Area": width * height * px * py, "Area_Unit": unit + '^2',
                        "Pixel_Width": px, "Pixel_Height": py, "Pixel_Unit": unit}
            return projection, metadata
        finally:
            self.release(selected)

    def save_projection(self, projection, path):
        """Reopen the actual TIFF and verify float values, dimensions and calibration."""
        if not self.FileSaver(projection).saveAsTiff(str(path)) or not path.is_file():
            raise RuntimeError(f"Failed to save 32-bit projection: {path}")
        reopened = self.IJ.openImage(str(path))
        try:
            if (reopened is None or int(reopened.getBitDepth()) != 32
                    or [int(v) for v in reopened.getDimensions()] != [int(v) for v in projection.getDimensions()]):
                raise RuntimeError(f"The stored projection is not the expected single-plane 32-bit TIFF: {path}")
            expected = self.np.asarray(projection.getProcessor().getPixels(), dtype=self.np.float32)
            actual = self.np.asarray(reopened.getProcessor().getPixels(), dtype=self.np.float32)
            if not self.np.array_equal(expected, actual):
                raise RuntimeError(f"TIFF export changed numeric projection intensities: {path}")
            a, b = projection.getCalibration(), reopened.getCalibration()
            if (not math.isclose(float(a.pixelWidth), float(b.pixelWidth), rel_tol=1e-6)
                    or not math.isclose(float(a.pixelHeight), float(b.pixelHeight), rel_tol=1e-6)
                    or str(a.getUnit()) != str(b.getUnit())):
                raise RuntimeError(f"TIFF export changed spatial calibration: {path}")
        finally:
            self.release(reopened)
        return sha256_file(path)

    def measure(self, projection, mask_path, limits):
        """ImageJ thresholds float pixels directly; no 8-bit intensity conversion."""
        processor = projection.getProcessor()
        processor.setThreshold(limits['lower'], limits['upper'], self.ImageProcessor.NO_LUT_UPDATE)
        values = self.np.asarray(processor.getPixels(), dtype=self.np.float32)
        expected = (values >= limits['lower']) & (values <= limits['upper'])
        mask_processor = processor.createMask()
        if mask_processor is None:
            raise RuntimeError("ImageJ could not create a mask from the explicit float threshold.")
        mask = self.ImagePlus(mask_path.name, mask_processor)
        mask.setCalibration(projection.getCalibration().copy())
        try:
            mask_values = self.np.asarray(mask_processor.getPixels(), dtype=self.np.int8).view(self.np.uint8)
            if not self.np.array_equal(mask_values, expected.astype(self.np.uint8) * 255):
                raise RuntimeError("ImageJ float-threshold mask disagrees with the raw projection values.")
            mask_processor.setThreshold(255.0, 255.0, self.ImageProcessor.NO_LUT_UPDATE)
            flags = int(self.Measurements.AREA) | int(self.Measurements.LIMIT) | int(self.Measurements.AREA_FRACTION)
            stats = mask.getStatistics(flags)
            positive = int(stats.pixelCount)
            total = int(values.size)
            percentage = float(stats.areaFraction)
            area = float(stats.area)
            cal = projection.getCalibration()
            if (positive != int(expected.sum())
                    or not math.isclose(percentage, positive / total * 100, abs_tol=1e-10, rel_tol=1e-10)
                    or not math.isclose(area, positive * float(cal.pixelWidth) * float(cal.pixelHeight), abs_tol=1e-9, rel_tol=1e-9)):
                raise RuntimeError("ImageJ area measurements failed reconciliation.")
            mask_processor.resetThreshold()
            if not self.FileSaver(mask).saveAsTiff(str(mask_path)):
                raise RuntimeError(f"Failed to save binary mask: {mask_path}")
            reopened = self.IJ.openImage(str(mask_path))
            try:
                if (reopened is None or int(reopened.getBitDepth()) != 8
                        or not self.np.array_equal(self.np.asarray(reopened.getProcessor().getPixels(), dtype=self.np.int8).view(self.np.uint8), mask_values)):
                    raise RuntimeError("Saved TIFF mask pixels differ from measured mask pixels.")
            finally:
                self.release(reopened)
            return {"FN_Positive_Pixels": positive, "FN_Area_Percent": percentage, "FN_Area": area,
                    "Threshold_Lower": limits['lower'],
                    "Threshold_Upper": None if limits['upper_unbounded'] else limits['upper'],
                    "Effective_Threshold_Upper": limits['upper'], "Threshold_Units": "Raw projection intensity",
                    "Requested_Threshold_Lower": limits['requested_lower'],
                    "Requested_Threshold_Upper": limits['requested_upper'], "Mask_Bit_Depth": 8,
                    "Mask_File_Name": mask_path.name, "Mask_Path": str(mask_path), "Mask_SHA256": sha256_file(mask_path)}
        finally:
            processor.resetThreshold()
            self.release(mask)

    def close(self):
        if self.ij is not None:
            try:
                self.ij.dispose()
            finally:
                self.ij = None


def verify_table(path, expected, fieldnames):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        saved = list(reader)
        if reader.fieldnames != fieldnames or len(saved) != len(expected):
            raise RuntimeError(f"CSV column or row count changed: {path.name}")
    for record, actual in zip(expected, saved):
        for key in fieldnames:
            value = record[key]
            if actual[key] != ('' if value is None else str(value)):
                raise RuntimeError(f"CSV export changed {key} for {record['Image_ID']}.")


def process_folder(source_folder, channel, method, limits, engine_holder, input_json, outputs):
    if not source_folder.is_dir():
        raise ValidationError(f"Source folder does not exist: {source_folder}")
    run_id, output = new_output_folder(source_folder)
    outputs.append(output)
    log = RunLog(output)
    mode = 'AREA_MEASUREMENT' if limits is not None else 'PROJECTIONS_ONLY'
    status = {"run_id": run_id, "script_version": SCRIPT_VERSION, "status": "RUNNING", "mode": mode,
              "started_utc": utc_now(), "source_folder": str(source_folder),
              "run_directory": str(output), "stage": "Original inventory", "processed_images": 0, "skipped_images": 0,
              "package_version": package_version()}
    parameters = {**status, "channel_index": channel, "projection_method": method.upper(), "projection_bit_depth": 32,
                  "intensity_scaling": "None", "resizing": "None; native XY resolution",
                  "background_subtraction": "None", "denoising": "None", "threshold": limits,
                  "threshold_policy": "Explicit inclusive raw float32 limits; no auto-threshold and no 8-bit rescaling",
                  "mask_values": [0, 255], "mask_bit_depth": 8, "denominator": "Full native-resolution XY image",
                  "selection_rule": "Every visible ND2/TIF/TIFF file directly in the source folder",
                  "input_json": str(input_json) if input_json else None,
                  "python": platform.python_version(), "platform": platform.platform(), "python_executable": sys.executable,
                  "imagej_endpoint": FIJI_ENDPOINT, "imagej_mode": "headless"}
    current_file = ''
    final_projection_manifest = final_summary = None
    try:
        save_json(output / 'run_status.json', status)
        save_json(output / 'run_parameters.json', parameters)
        log.event('STARTED', 'Run', f"Version {SCRIPT_VERSION}; {mode}; {method.upper()}32; channel {channel}")
        log.event('INFO', 'Rules', 'Original intensities and XY resolution are preserved. SUM also sums background. No denoising is applied.')
        log.event('INFO', 'Threshold', json.dumps(limits) if limits else 'Not supplied: projections only; no FN area table or masks will be created.')
        status['stage'] = 'Original inventory'
        selected = original_inventory(source_folder)
        save_csv(output / 'input_manifest.csv', MANIFEST_COLUMNS, selected)
        log.event('INFO', 'Inventory', f"{len(selected)} original image(s); subfolders and hidden files are not processed.")
        for index, entry in enumerate(selected, 1):
            log.event('INFO', 'Fingerprint', f"{index}/{len(selected)} {entry['File_Name']}")
            entry['SHA256'] = sha256_file(Path(entry['Path']))
        save_csv(output / 'input_manifest.csv', MANIFEST_COLUMNS, selected)
        status.update(input_images=len(selected), stage='ImageJ initialization')
        save_json(output / 'run_status.json', status)
        if engine_holder[0] is None:
            engine_holder[0] = ImageJEngine(log)
        engine = engine_holder[0]
        parameters['runtime_versions'] = engine.versions
        save_json(output / 'run_parameters.json', parameters)
        projections = output / 'Projections_32bit'
        projections.mkdir()
        masks = output / 'Masks'
        if limits is not None:
            masks.mkdir()
        projection_partial = output / 'FN_Projection_Manifest.partial.csv'
        summary_partial = output / 'Fibronectin_Area_Summary.partial.csv'
        projection_rows, summary_rows = [], []
        for index, entry in enumerate(selected, 1):
            path = Path(entry['Path'])
            current_file = path.name
            status.update(stage='Original projection', current_file=current_file)
            save_json(output / 'run_status.json', status)
            log.event('INFO', 'Image', f"{index}/{len(selected)} {path.name}")
            original = projection = None
            try:
                original, reader_name = engine.open_original(path, channel)
                projection, measured = engine.project(original, channel, method)
                engine.release(original)
                original = None
                projection_path = projections / f"{path.name}_FN_{method.upper()}32.tif"
                projection_hash = engine.save_projection(projection, projection_path)
                row = {'File_Name': path.name, 'Image_ID': entry['Image_ID'], **measured,
                       'Source_Original_Path': str(path), 'Source_Projection_Path': str(projection_path),
                       'Projection_SHA256': projection_hash, 'Source_Folder': str(source_folder),
                       'Source_Reader': reader_name, 'Source_SHA256': entry['SHA256'], 'Program_Version': SCRIPT_VERSION, 'Run_ID': run_id}
                projection_rows.append(row)
                save_csv(projection_partial, PROJECTION_COLUMNS, projection_rows)
                if limits is not None:
                    mask_path = masks / ('FN_Mask_' + projection_path.name)
                    area = engine.measure(projection, mask_path, limits)
                    summary_rows.append({**row, **area})
                    save_csv(summary_partial, SUMMARY_COLUMNS, summary_rows)
                    log.event('PASS', 'Area', f"{area['FN_Positive_Pixels']}/{measured['Total_Pixels']} positive pixels ({area['FN_Area_Percent']:.6f}%).")
                log.event('PASS', 'Projection', f"{method.upper()}32; {measured['Number_of_Z_Stacks']} Z slices; "
                          f"raw range {measured['Projection_Min']:.6g} to {measured['Projection_Max']:.6g}; TIFF values verified.")
                status['processed_images'] = index
                save_json(output / 'run_status.json', status)
            finally:
                engine.release(projection)
                engine.release(original)
        status['stage'] = 'Final verification'
        current_selected = original_inventory(source_folder)
        if [entry['Path'] for entry in current_selected] != [entry['Path'] for entry in selected]:
            raise ValidationError('The original image inventory changed during processing.')
        for entry in selected:
            if sha256_file(Path(entry['Path'])) != entry['SHA256']:
                raise ValidationError('An original image changed during processing: ' + entry['File_Name'])
        z_counts = sorted({row['Number_of_Z_Stacks'] for row in projection_rows})
        if method == 'sum' and len(z_counts) > 1:
            log.event('WARNING', 'SUM comparability', f"Different Z counts: {z_counts}. A common SUM threshold is Z-count dependent. No automatic normalization was applied.")
        verify_table(projection_partial, projection_rows, PROJECTION_COLUMNS)
        if len(projection_rows) != len(selected):
            raise RuntimeError('Not every input image received an original projection.')
        if limits is not None:
            verify_table(summary_partial, summary_rows, SUMMARY_COLUMNS)
            if len(summary_rows) != len(selected):
                raise RuntimeError('Area summary does not include every selected image.')
        parameters.update(runtime_versions=engine.versions, z_slice_counts=z_counts, included_images=len(selected))
        save_json(output / 'run_parameters.json', parameters)
        final_projection_manifest = output / 'FN_Projection_Manifest.csv'
        projection_partial.rename(final_projection_manifest)
        if limits is not None:
            final_summary = output / 'Fibronectin_Area_Summary.csv'
            summary_partial.rename(final_summary)
        status.update(status='SUCCESS', stage='Completed', ended_utc=utc_now(),
                      generated_projections=len(projection_rows), generated_masks=len(summary_rows),
                      projection_manifest=str(final_projection_manifest), summary=str(final_summary) if final_summary else None,
                      failed_images=0, unprocessed_images=0)
        status.pop('current_file', None)
        save_json(output / 'run_status.json', status)
        log.event('SUCCESS', 'Run', f"{len(projection_rows)} original images projected; {len(summary_rows)} area measurements. Output: {output}")
        return True, output
    except (Exception, KeyboardInterrupt) as error:
        state = 'CANCELLED' if isinstance(error, KeyboardInterrupt) else 'VALIDATION_FAILED' if isinstance(error, ValidationError) else 'ERROR'
        message = str(error) or 'Run interrupted.'
        for published in (final_summary, final_projection_manifest):
            if published is not None and published.exists():
                published.rename(published.with_name(published.stem + '.partial.csv'))
        log.event('ERROR', status['stage'], message)
        save_csv(output / 'errors.csv', ['File_Name', 'Stage', 'Issue'],
                 [{'File_Name': current_file, 'Stage': status['stage'], 'Issue': message}])
        (output / 'traceback.txt').write_text(traceback.format_exc(), encoding='utf-8')
        failed_images = int(bool(current_file) and status['stage'] == 'Original projection')
        status.update(status=state, ended_utc=utc_now(), error=message, failed_images=failed_images,
                      unprocessed_images=max(0, status.get('input_images', 0) - status['processed_images'] - failed_images))
        save_json(output / 'run_status.json', status)
        log.event('FAILED', 'Run', f"Partial results and diagnostics retained: {output}")
        if isinstance(error, KeyboardInterrupt):
            raise
        return False, output
    finally:
        log.close()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Measure FN area from original ND2/TIFF images using native-resolution SUM32 projections.")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {package_version()} (area script {SCRIPT_VERSION})")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-i", "--input", help="JSON containing folder_paths; relative JSON paths use the working directory.")
    source.add_argument("--folder", help="One original-image folder; relative paths use the working directory.")
    parser.add_argument("--channel", type=int, help="Fibronectin channel, numbered from 1; prompt once if omitted.")
    parser.add_argument("--projection", choices=["sum", "mean", "max"], default="sum",
                        help="Default: sum. Every output is a native-resolution 32-bit TIFF.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("-t", "--threshold", nargs="+", metavar="VALUE",
                      help="Inclusive LOWER [UPPER] in raw projection units. Default LOWER: 2000. "
                           "Omit UPPER or use inf for the largest finite float32 value.")
    mode.add_argument("--projections-only", action="store_true",
                      help="Save 32-bit projections without masks or an area summary.")
    return parser.parse_args(argv)


def finish_imagej(holder, outputs):
    """Dispose the context and the same ImageJ worker pool as thickness."""
    errors = []
    try:
        if holder[0] is not None:
            holder[0].close()
    except (Exception, KeyboardInterrupt):
        errors.append(traceback.format_exc())
    try:
        if __package__:
            from .cli import _shutdown_imagej_workers
        else:
            from cli import _shutdown_imagej_workers
        _shutdown_imagej_workers()
    except (Exception, KeyboardInterrupt):
        errors.append(traceback.format_exc())

    for output in outputs:
        log = None
        try:
            log = RunLog(output, append=True)
            if errors:
                log.event("ERROR", "ImageJ shutdown", "Context or worker shutdown failed; see shutdown_traceback.txt.")
                (output / "shutdown_traceback.txt").write_text("\n".join(errors), encoding="utf-8")
                status_path = output / "run_status.json"
                if status_path.exists():
                    status = json.loads(status_path.read_text(encoding="utf-8"))
                    status.update(shutdown_status="ERROR", ended_utc=utc_now())
                    if status["status"] == "SUCCESS":
                        status.update(status="ERROR", stage="ImageJ shutdown")
                    save_json(status_path, status)
            else:
                log.event("INFO", "ImageJ shutdown", "ImageJ context and workers closed; returning to the terminal.")
        except OSError as error:
            errors.append(f"Could not record shutdown in {output}: {error}")
        finally:
            if log is not None:
                log.close()
    if errors:
        print("ImageJ shutdown failed:\n" + "\n".join(errors), file=sys.stderr, flush=True)
    return not errors


def main(argv=None):
    args = parse_args(argv)
    holder, outputs, folders = [None], [], []
    exit_code = 0
    try:
        folders, input_json = read_source_folders(args)
        bounds = args.threshold if args.threshold is not None else [DEFAULT_THRESHOLD_LOWER]
        if not 1 <= len(bounds) <= 2:
            raise ValidationError("Use -t LOWER or -t LOWER UPPER.")
        lower = bounds[0]
        upper = bounds[1] if len(bounds) == 2 else None
        limits = None if args.projections_only else threshold_settings(lower, upper)
        channel = args.channel
        if channel is None:
            try:
                channel = int(input("Enter fibronectin channel index (starting from 1): ").strip())
            except EOFError as error:
                raise ValidationError("Supply --channel when interactive input is unavailable.") from error
        if isinstance(channel, bool) or not isinstance(channel, int) or channel < 1:
            raise ValidationError("The fibronectin channel must be an integer greater than or equal to 1.")
        print(f"UMA-tools {package_version()}; area {SCRIPT_VERSION}: {args.projection.upper()}32; channel {channel}; "
              + (f"threshold {limits['lower']:.9g} to {limits['upper']:.9g}."
                 if limits is not None else "projections only."), flush=True)
        completed = failed = 0
        for folder in folders:
            try:
                ok, _ = process_folder(folder, channel, args.projection, limits, holder, input_json, outputs)
                completed += int(ok)
                failed += int(not ok)
            except (ValidationError, OSError) as error:
                failed += 1
                print(f"Cannot process {folder}: {error}", file=sys.stderr, flush=True)
                output = save_startup_error(error, [folder], args)
                if output is not None:
                    outputs.append(output)
        print(f"Finished: {completed} successful folder(s), {failed} failed folder(s).", flush=True)
        exit_code = 0 if failed == 0 else 1
    except KeyboardInterrupt as error:
        print("Processing cancelled.", file=sys.stderr, flush=True)
        if not outputs:
            output = save_startup_error(error, folders, args)
            if output is not None:
                outputs.append(output)
        exit_code = 130
    except Exception as error:
        print(f"Cannot start: {error}", file=sys.stderr, flush=True)
        output = save_startup_error(error, folders, args)
        if output is not None:
            outputs.append(output)
        exit_code = 2 if isinstance(error, (ValidationError, OSError, ValueError)) else 1
    finally:
        if not finish_imagej(holder, outputs) and exit_code == 0:
            exit_code = 1
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
