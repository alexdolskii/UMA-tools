#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
from datetime import datetime
from pathlib import Path

import pandas as pd
from scyjava import jimport

from .config import load_json
from .imagej import (
    initialize_imagej,
)
from .run import scoped_file_log, unique_output

_LOGGER = logging.getLogger(__name__)


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
    """
    Prompt the user to choose the file type (.nd2 or .tiff).

    Returns:
        str: The file extension chosen by the user.
    """
    print("\nSelect the file type to process:")
    print("1. .nd2")
    print("2. .tiff")
    choice = input("Enter 1 for .nd2 or 2 for .tiff: ").strip()
    if choice == "1":
        return ".nd2"
    elif choice == "2":
        return ".tiff"
    else:
        raise ValueError(
            "Invalid choice. Please run the script again and select 1 or 2."
        )


def get_fibronectin_channel():
    """
    Prompt for the 1-based fibronectin channel index.

    Channel availability is checked against each image during
    processing.

    Returns:
        int: The fibronectin channel index.
    """
    try:
        val = int(
            input(
                "Enter fibronectin channel index (starting from 1): "
            ).strip()
        )
    except ValueError:
        raise ValueError(
            "Please enter an integer for the fibronectin channel index."
        )
    if val < 1:
        raise ValueError("The fibronectin channel index must be at least 1.")
    return val


def get_folder_paths(input_file_path):
    """
    Reads an input JSON file containing folder paths.

    Args:
        input_file_path (str): Path to the input JSON file.

    Returns:
        List[str]: List of valid folder paths.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file does not contain folder paths or no
        valid
            folders.
    """
    if os.path.basename(input_file_path).startswith("._"):
        raise ValueError(
            f"macOS metadata files cannot be used as input: {input_file_path}"
        )
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"File '{input_file_path}' does not exist.")

    if not input_file_path.lower().endswith(".json"):
        raise ValueError("Input file must be in .json format.")

    data = load_json(
        Path(input_file_path), encoding="utf-8", reject_metadata=False
    )

    folder_paths = data.get("folder_paths", [])
    if not folder_paths:
        raise ValueError("Input file does not contain folder paths.")

    valid_folder_paths = []
    for folder_path in folder_paths:
        if os.path.isdir(folder_path):
            files = [
                f
                for f in os.listdir(folder_path)
                if not f.startswith(".")
                and os.path.isfile(os.path.join(folder_path, f))
            ]
            num_files = len(files)
            file_types = set(
                [
                    os.path.splitext(f)[1].lower()
                    for f in files
                    if not f.startswith(".")
                ]
            )
            print(f"\nFolder: {folder_path}")
            print(f"Number of files: {num_files}")
            print(f"File types: {', '.join(file_types)}")
            valid_folder_paths.append(folder_path)
        else:
            print(f"\nFolder '{folder_path}' does not exist.")

    if not valid_folder_paths:
        raise ValueError("No available folders for processing.")

    print(
        f"\nFound {len(valid_folder_paths)} available folders for processing."
    )
    return valid_folder_paths


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
    print("  Opening image...")
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
    print("  Extracting fibronectin channel...")
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
    print("  Performing Reslice...")
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
    print("  Performing Z projection...")
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
    print("  Applying Maximum filter...")
    IJ.run(projected_imp, "Maximum...", "radius=2")

    print("  Applying Gaussian Blur...")
    IJ.run(projected_imp, "Gaussian Blur...", "sigma=2 scaled")

    print("  Subtracting background...")
    IJ.run(projected_imp, "Subtract Background...", "rolling=50 sliding")

    print("  Applying threshold...")
    IJ.setAutoThreshold(projected_imp, "Otsu dark no-reset")
    IJ.run("Options...", "black")
    IJ.run(projected_imp, "Convert to Mask", "")


def calculate_local_thickness(IJ, projected_imp, filename: str):
    """Run the masked, calibrated, silent Local Thickness plugin."""
    # Run Local Thickness
    print("  Running Local Thickness...")
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
    print("  Measuring thickness...")
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

    # Extract results
    if rt is None or rt.getCounter() == 0:
        _LOGGER.warning("No measurements.")
        area = None
        std_dev = None
        min_thickness = None
        max_thickness = None
        median_thickness = None
    else:
        try:
            row = rt.getCounter() - 1
            area = rt.getValue("Area", row)
            std_dev = rt.getValue("StdDev", row)
            min_thickness = rt.getValue("Min", row)
            max_thickness = rt.getValue("Max", row)
            median_thickness = rt.getValue("Median", row)
            _LOGGER.info(
                f"  Results - Area: {area}, StdDev: {std_dev}, "
                f"Min: {min_thickness}, Max: {max_thickness}, "
                f"Median: {median_thickness}"
            )
        except Exception as e:
            _LOGGER.error(f"Error reading measurements: {e}")
            area = None
            std_dev = None
            min_thickness = None
            max_thickness = None
            median_thickness = None
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
    print("  Closed all images.\n")

    return {"File_Name": filename, **measurements}


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
    image_files = [
        f
        for f in os.listdir(folder)
        if f.lower().endswith(file_extension)
        and not f.startswith(".")
        and os.path.isfile(os.path.join(folder, f))
    ]

    if len(image_files) == 0:
        print(
            f"No '{file_extension}' files found in folder '{folder}'. "
            "Skipping."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    _, output = unique_output(
        Path(folder),
        "Thickness_assay_results_",
        timestamp=timestamp,
        include_pid=True,
    )
    results_folder = str(output)

    with scoped_file_log(_LOGGER, output):
        try:
            _LOGGER.info(f"\nProcessing folder: {folder}")
            _LOGGER.info(f"Results will be saved in: {results_folder}")

            summary_data = []
            for filename in image_files:
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
                if result is not None:
                    summary_data.append(result)

            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_file_path = os.path.join(
                    results_folder, "Thickness_Summary.csv"
                )
                summary_df.to_csv(summary_file_path, index=False)
                _LOGGER.info(
                    f"Folder analysis complete. "
                    f"Data saved to '{summary_file_path}'."
                )
            else:
                _LOGGER.warning(
                    f"No data to save in summary for folder '{folder}'."
                )
        except BaseException:
            _LOGGER.exception("Thickness folder analysis failed.")
            raise


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
    for folder in folder_paths:
        # Run analysis
        process_single_folder(
            ij,
            IJ,
            WindowManager,
            Duplicator,
            ResultsTable,
            folder,
            file_extension,
            fibronectin_channel,
        )
    print("\nAll folders have been processed.")


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
    # Setting up logging
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    _LOGGER.info("Thickness implementation: %s", Path(__file__).resolve())

    folder_paths = get_folder_paths(input_json_path)
    ij = initialize_imagej()
    try:
        (IJ, ImagePlus, WindowManager, ResultsTable, Duplicator, System) = (
            import_java_classes()
        )
        file_extension = get_file_type_choice()
        fibronectin_channel = get_fibronectin_channel()

        start_analysis = (
            input("\nDo you want to start processing? (y/n): ").strip().lower()
        )
        if start_analysis in ("no", "n"):
            raise ValueError("Analysis canceled by user.")
        elif start_analysis not in ("yes", "y", "no", "n"):
            raise ValueError("Incorrect input. Please enter y/n or yes/no")

        process_all_folders(
            ij,
            IJ,
            WindowManager,
            Duplicator,
            ResultsTable,
            folder_paths,
            file_extension,
            fibronectin_channel,
        )
    except Exception:
        _LOGGER.exception("Thickness analysis failed.")
        raise
    finally:
        print("Disposing of ImageJ context...")
        ij.dispose()
    print("Thickness analysis is successfully completed.")
