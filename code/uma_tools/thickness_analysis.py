#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import math
import os
from datetime import datetime
from pathlib import Path

import pandas as pd
from scyjava import jimport

from .config import read_config
from .files import assay_directory
from .image_run import ImageRun, image_names
from .imagej import (
    initialize_imagej,
)
from .progress import (
    CANCELLATIONS,
    ask_channel,
    ask_choice,
    confirm_start,
    folder_error,
    folder_logged,
    folder_scope,
    outcome,
    phase,
    register_sources,
)
from .run import scoped_file_log, unique_output
from .runtime import temporary_path

_LOGGER = logging.getLogger(__name__)

PROGRESS_STAGES = ("Thickness",)
PROGRESS_OPERATIONS = (
    "Opening image...",
    "Extracting fibronectin channel...",
    "Performing Reslice...",
    "Performing Z projection...",
    "Applying Maximum filter...",
    "Applying Gaussian Blur...",
    "Subtracting background...",
    "Applying threshold...",
    "Running Local Thickness...",
    "Measuring thickness...",
    "Closed all images.",
)


def import_java_classes():
    """
    Import necessary Java classes for image processing.

    Returns:
        tuple: A tuple containing references to imported classes.
    """
    IJ = jimport("ij.IJ")
    ImagePlus = jimport("ij.ImagePlus")
    WindowManager = jimport("ij.WindowManager")
    ResultsTable = jimport("ij.measure.ResultsTable")
    Duplicator = jimport("ij.plugin.Duplicator")
    System = jimport("java.lang.System")
    return IJ, ImagePlus, WindowManager, ResultsTable, Duplicator, System


def get_file_type_choice():
    """Retry unsupported choices without restarting the command."""
    return ask_choice(
        "Enter 1 for .nd2 or 2 for .tiff: ", {"1": ".nd2", "2": ".tiff"}
    )


def get_fibronectin_channel():
    return ask_channel()


def get_folder_paths(input_file_path):
    return [str(path) for path in read_config(Path(input_file_path))]


def open_fibronectin_channel(
    ij,
    IJ,
    Duplicator,
    file_path: str,
    filename: str,
    fibronectin_channel: int,
):
    """
    Open the original image and duplicate the requested full Z/T
    channel.
    """
    # Open image
    phase("Opening image...")
    img = ij.io().open(file_path)
    if img is None:
        _LOGGER.warning(f"Failed to open image: {filename}")
        IJ.run("Close All")
        return

    imp = ij.convert().convert(img, jimport("ij.ImagePlus"))
    if imp is None:
        _LOGGER.warning(f"Failed to convert image '{filename}' to ImagePlus.")
        IJ.run("Close All")
        return

    # Extract fibronectin channel
    phase("Extracting fibronectin channel...")
    if fibronectin_channel > imp.getNChannels() or fibronectin_channel < 1:
        _LOGGER.warning(f"Invalid fibronectin channel for image {filename}")
        imp.close()
        IJ.run("Close All")
        return

    imp_fibronectin = Duplicator().run(
        imp,
        fibronectin_channel,
        fibronectin_channel,
        1,
        imp.getNSlices(),
        1,
        imp.getNFrames(),
    )
    imp.close()

    if imp_fibronectin is None:
        _LOGGER.warning(
            f"Failed to extract fibronectin channel in {filename}."
        )
        IJ.run("Close All")
        return

    imp_fibronectin.setTitle(f"{filename}_C{fibronectin_channel}")
    return imp_fibronectin


def reslice_and_project(IJ, imp_fibronectin, filename: str):
    """
    Apply the original reslice macro and calibrated MAX Z projection.
    """
    # Reslice
    phase("Performing Reslice...")
    # Batch mode preserves the original macro options without GUI
    # windows.
    interpreter = jimport("ij.macro.Interpreter")()
    resliced_imp = interpreter.runBatchMacro(
        'run("Reslice [/]...", "output=0.500 start=Top flip rotate avoid");',
        imp_fibronectin,
    )
    imp_fibronectin.close()
    if resliced_imp is None:
        _LOGGER.warning(f"Failed to perform Reslice for {filename}")
        IJ.run("Close All")
        return
    resliced_imp.setTitle(f"Reslice_of_{filename}")
    calibration = resliced_imp.getCalibration()
    _LOGGER.info(
        "Reslice calibration: pixelWidth=%s, pixelHeight=%s, unit=%s",
        calibration.pixelWidth,
        calibration.pixelHeight,
        calibration.getUnit(),
    )

    # Z Project
    phase("Performing Z projection...")
    # Match the original Z Project macro: preserve spatial calibration
    # and
    # project the current time frame. The low-level doProjection() loses
    # both.
    projected_imp = jimport("ij.plugin.ZProjector").run(resliced_imp, "max")
    resliced_imp.close()
    if projected_imp is None:
        _LOGGER.warning(f"Failed to perform Z projection for {filename}")
        IJ.run("Close All")
        return
    projected_imp.setTitle(f"MAX_Reslice_of_{filename}")
    calibration = projected_imp.getCalibration()
    _LOGGER.info(
        "Projection calibration: pixelWidth=%s, pixelHeight=%s, unit=%s",
        calibration.pixelWidth,
        calibration.pixelHeight,
        calibration.getUnit(),
    )
    return projected_imp


def create_thickness_mask(IJ, projected_imp) -> None:
    """
    Run the original filters and Otsu mask in their unchanged order.
    """
    # Filters and threshold
    phase("Applying Maximum filter...")
    IJ.run(projected_imp, "Maximum...", "radius=2")

    phase("Applying Gaussian Blur...")
    IJ.run(projected_imp, "Gaussian Blur...", "sigma=2 scaled")

    phase("Subtracting background...")
    IJ.run(projected_imp, "Subtract Background...", "rolling=50 sliding")

    phase("Applying threshold...")
    IJ.setAutoThreshold(projected_imp, "Otsu dark no-reset")
    IJ.run("Options...", "black")
    IJ.run(projected_imp, "Convert to Mask", "")


def calculate_local_thickness(IJ, projected_imp, filename: str):
    """Run the masked, calibrated, silent Local Thickness plugin."""
    # Run Local Thickness
    phase("Running Local Thickness...")
    # This is the same plugin registered by the masked/calibrated/silent
    # menu
    # command, but processImage returns its result without opening a
    # window.
    local_thickness = jimport("sc.fiji.localThickness.LocalThicknessWrapper")()
    local_thickness.setSilence(True)
    local_thickness.setShowOptions(False)
    local_thickness.maskThicknessMap = True
    local_thickness.calibratePixels = True
    local_thickness_imp = local_thickness.processImage(projected_imp)

    if local_thickness_imp is None:
        _LOGGER.warning("  Could not retrieve Local Thickness image.")
        IJ.run("Close All")
        return

    local_thickness_imp.setTitle(f"Local_Thickness_of_{filename}")
    IJ.run(local_thickness_imp, "Fire", "")
    return local_thickness_imp


def measure_thickness(ResultsTable, local_thickness_imp) -> dict:
    """Measure five statistics using an image-private results table."""
    # Measure the explicit image into a private table without a Results
    # window.
    # Preserve the original area, standard deviation, min/max and median
    # flags.
    phase("Measuring thickness...")
    measurements = jimport("ij.measure.Measurements")
    flags = (
        measurements.AREA
        | measurements.STD_DEV
        | measurements.MIN_MAX
        | measurements.MEDIAN
    )
    rt = ResultsTable()
    rt.setPrecision(3)
    analyzer = jimport("ij.plugin.filter.Analyzer")(
        local_thickness_imp, flags, rt
    )
    analyzer.measure()

    if rt is None or rt.getCounter() == 0:
        raise ValueError("No thickness measurements returned")
    row = rt.getCounter() - 1
    area = rt.getValue("Area", row)
    std_dev = rt.getValue("StdDev", row)
    min_thickness = rt.getValue("Min", row)
    max_thickness = rt.getValue("Max", row)
    median_thickness = rt.getValue("Median", row)
    if not all(
        math.isfinite(value)
        for value in (
            area,
            std_dev,
            min_thickness,
            max_thickness,
            median_thickness,
        )
    ):
        raise ValueError("Thickness measurements contain nonfinite values")
    _LOGGER.info(
        "Area=%s, StdDev=%s, Min=%s, Max=%s, Median=%s",
        area,
        std_dev,
        min_thickness,
        max_thickness,
        median_thickness,
    )
    return {
        "Area": area,
        "StdDev": std_dev,
        "Min": min_thickness,
        "Max": max_thickness,
        "Median": median_thickness,
    }


def process_single_file(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder,
    filename,
    fibronectin_channel,
    results_folder,
) -> dict:
    """
    Process a single image file.

    Returns:
        dict: A dictionary with measurement results for this file.
    """
    file_path = os.path.join(folder, filename)
    _LOGGER.info(f"Processing file: {filename}")

    imp_fibronectin = projected_imp = local_thickness_imp = None
    try:
        imp_fibronectin = open_fibronectin_channel(
            ij, IJ, Duplicator, file_path, filename, fibronectin_channel
        )
        if imp_fibronectin is None:
            return None
        projected_imp = reslice_and_project(IJ, imp_fibronectin, filename)
        if projected_imp is None:
            return None
        create_thickness_mask(IJ, projected_imp)

        # Save mask image
        mask_image_path = os.path.join(results_folder, f"Mask_{filename}.tif")
        IJ.saveAs(projected_imp, "Tiff", mask_image_path)
        _LOGGER.info(f"Mask saved to '{mask_image_path}'.")

        local_thickness_imp = calculate_local_thickness(
            IJ, projected_imp, filename
        )
        if local_thickness_imp is None:
            return None
        measurements = measure_thickness(ResultsTable, local_thickness_imp)

        # Save thickness image
        thickness_path = os.path.join(
            results_folder, f"Local_Thickness_{filename}.tif"
        )
        IJ.saveAs(local_thickness_imp, "Tiff", thickness_path)
        _LOGGER.info(f"Local Thickness image saved to '{thickness_path}'.")

        projected_imp.close()
        local_thickness_imp.close()
        IJ.run("Close All")
        phase("Closed all images.\n")

        return {"File_Name": filename, **measurements}
    finally:
        for image in (local_thickness_imp, projected_imp, imp_fibronectin):
            if image is not None:
                image.close()
        IJ.run("Close All")


@folder_logged("folder")
def process_single_folder(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder,
    file_extension,
    fibronectin_channel,
):
    """
    Process all image files in a single folder.

    Args:
        ij (imagej.ImageJ): The ImageJ instance.
        IJ, WindowManager, Duplicator, ResultsTable: Java class
        references.
        folder (str): The folder path to process.
        file_extension (str): The file extension to process (.nd2 or
        .tiff).
        fibronectin_channel (int): The fibronectin channel index.
    """
    image_files = image_names(folder, (file_extension,))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    _, output = unique_output(
        assay_directory(Path(folder), create=True),
        "Thickness_assay_results_",
        timestamp=timestamp,
        include_pid=True,
    )
    results_folder = str(output)

    run = ImageRun(
        output,
        folder,
        image_files,
        {"channel_index": fibronectin_channel, "extension": file_extension},
        stages=PROGRESS_STAGES,
        operations=PROGRESS_OPERATIONS,
    )
    with scoped_file_log(_LOGGER, output), run:
        summary_data = []
        summary_path = output / "Thickness_Summary.csv"
        _LOGGER.info("Results: %s", output)
        for filename in image_files:
            with run.attempt(filename, "Thickness", final=True):
                result = process_single_file(
                    ij,
                    IJ,
                    WindowManager,
                    Duplicator,
                    ResultsTable,
                    folder,
                    filename,
                    fibronectin_channel,
                    results_folder,
                )
                if result is None:
                    raise ValueError(
                        "Image processing returned no result; see log.log "
                        "for the failed operation"
                    )
                if not all(
                    math.isfinite(float(result[key]))
                    for key in ("Area", "StdDev", "Min", "Max", "Median")
                ):
                    raise ValueError("Invalid thickness measurements")
                # Commit a row only after its summary has been saved.
                candidate = summary_data + [result]
                pending = temporary_path(
                    summary_path, summary_path.with_suffix(".pending.csv")
                )
                pd.DataFrame(candidate).to_csv(pending, index=False)
                pending.replace(summary_path)
                summary_data = candidate
        return run.finish(summary_path)


def process_all_folders(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder_paths,
    file_extension,
    fibronectin_channel,
):
    """
    Process all provided folders.

    Args:
        ij (imagej.ImageJ): The ImageJ instance.
        IJ, WindowManager, Duplicator, ResultsTable: Java class
        references.
        folder_paths (list[str]): List of folder paths to process.
        file_extension (str): The file extension to process (.nd2 or
        .tiff).
        fibronectin_channel (int): The fibronectin channel index.
    """
    failed = 0
    for folder in folder_paths:
        try:
            state = process_single_folder(
                ij,
                IJ,
                WindowManager,
                Duplicator,
                ResultsTable,
                folder,
                file_extension,
                fibronectin_channel,
            )
            failed += state != "SUCCESS"
        except CANCELLATIONS:
            raise
        except Exception as error:
            _LOGGER.exception("Folder failed: %s", folder)
            outcome("FAILED", f"{folder}: {error}")
            failed += 1
    return failed


def main(input_json_path: str) -> None:
    """
    The main function to perform thickness analysis on
    image files stored in specified folders. It initializes ImageJ,
    reads input parameters from a JSON file (specified
    via command-line arguments), prompts the user for file type and
    channel information, processes images in each folder,
    applies various filters and measurements,
    and saves the resulting data and images.

    Args:
        input_json_path: path to a json file
    """
    folder_paths = get_folder_paths(input_json_path)
    register_sources(folder_paths, input_json_path)
    file_extension = get_file_type_choice()
    fibronectin_channel = get_fibronectin_channel()
    confirm_start()
    ij = None
    failed = 0
    try:
        for folder in folder_paths:
            try:
                with folder_scope(folder):
                    classes = (None,) * 6
                    if image_names(folder, (file_extension,)):
                        if ij is None:
                            phase("Starting Fiji")
                            ij = initialize_imagej()
                        classes = import_java_classes()
                    IJ, _, WindowManager, ResultsTable, Duplicator, _ = classes
                    state = process_single_folder(
                        ij,
                        IJ,
                        WindowManager,
                        Duplicator,
                        ResultsTable,
                        folder,
                        file_extension,
                        fibronectin_channel,
                    )
                    failed += state != "SUCCESS"
            except CANCELLATIONS:
                raise
            except Exception as error:
                folder_error(folder, error)
                failed += 1
    finally:
        if ij is not None:
            phase("Closing ImageJ context")
            ij.dispose()
    outcome(
        "FINISHED",
        f"Thickness: {len(folder_paths) - failed} complete; "
        f"{failed} incomplete or failed folder(s).",
    )
    return int(bool(failed))
