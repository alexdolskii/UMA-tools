#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script performs thickness analysis on image files stored in specified
folders. It initializes ImageJ, reads input parameters from a JSON file
(specified via command-line arguments), prompts the user for file type and
channel information, processes images in each folder, applies various filters
and measurements, and saves the resulting data and images.
"""

import argparse
import imagej
import json
import os
import pandas as pd
import sys
import time
from datetime import datetime
from pathlib import Path
from scyjava import jimport


def initialize_imagej():
    """
    Initialize ImageJ in headless mode.

    Returns:
        ij (imagej.ImageJ): The initialized ImageJ instance.
    """
    # Attempt to initialize ImageJ headless mode
    ij = imagej.init('sc.fiji:fiji', mode='headless')
    return ij


def import_java_classes():
    """
    Import necessary Java classes for image processing.

    Returns:
        tuple: A tuple containing references to imported classes.
    """
    IJ = jimport('ij.IJ')
    ImagePlus = jimport('ij.ImagePlus')
    WindowManager = jimport('ij.WindowManager')
    ResultsTable = jimport('ij.measure.ResultsTable')
    Duplicator = jimport('ij.plugin.Duplicator')
    System = jimport('java.lang.System')
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
    if choice == '1':
        return '.nd2'
    elif choice == '2':
        return '.tiff'
    else:
        raise ValueError(
            "Invalid choice. Please run the script again and select 1 or 2."
        )


def get_num_channels():
    """
    Prompt the user for the number of channels in the files.

    Returns:
        int: The number of channels.
    """
    while True:
        val = input("Enter the number of channels in the files: ").strip()
        if val.isdigit():
            num = int(val)
            if num < 1:
                print("The number of channels must be at least 1.")
            else:
                return num
        else:
            print("Please enter an integer for the number of channels.")


def get_fibronectin_channel(num_channels):
    """
    Prompt the user for the fibronectin channel number if there are multiple
    channels.

    Args:
        num_channels (int): The total number of channels.

    Returns:
        int: The fibronectin channel number.
    """
    if num_channels == 1:
        return 1
    while True:
        val = input(
            f"Enter the channel number (1-{num_channels}) representing "
            "fibronectin: "
        ).strip()
        if val.isdigit():
            fib_ch = int(val)
            if 1 <= fib_ch <= num_channels:
                return fib_ch
            else:
                print(
                    f"Channel number must be between 1 and {num_channels}."
                )
        else:
            print("Please enter an integer for the channel number.")


def get_folder_paths(input_file_path):
    """
    Reads an input JSON file containing folder paths.

    Args:
        input_file_path (str): Path to the input JSON file.

    Returns:
        List[str]: List of valid folder paths.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file does not contain folder paths or no valid
            folders.
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"File '{input_file_path}' does not exist.")

    if not input_file_path.lower().endswith('.json'):
        raise ValueError("Input file must be in .json format.")

    with open(input_file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    folder_paths = data.get('folder_paths', [])
    if not folder_paths:
        raise ValueError("Input file does not contain folder paths.")

    valid_folder_paths = []
    for folder_path in folder_paths:
        if os.path.isdir(folder_path):
            files = os.listdir(folder_path)
            num_files = len(files)
            file_types = set([
                os.path.splitext(f)[1].lower()
                for f in files if not f.startswith(".")
                ])
            print(f"\nFolder: {folder_path}")
            print(f"Number of files: {num_files}")
            print(f"File types: {', '.join(file_types)}")
            valid_folder_paths.append(folder_path)
        else:
            print(f"\nFolder '{folder_path}' does not exist.")

    if not valid_folder_paths:
        raise ValueError("No available folders for processing.")

    print(f"\nFound {len(valid_folder_paths)} available folders "
          "for processing.")
    return valid_folder_paths


def confirm_processing():
    """
    Ask the user if they want to proceed with processing.

    Returns:
        bool: True if user wants to proceed, False otherwise.
    """
    proceed = input("\nDo you want to start processing? (y/n): ").strip().lower()
    return proceed == 'y'

def process_single_file(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder,
    filename,
    fibronectin_channel,
    file_extension,
    results_folder
):
    """
    Process a single image file.

    Returns:
        dict: A dictionary with measurement results for this file.
    """
    file_path = os.path.join(folder, filename)
    print(f"Processing file: {filename}")

    # Open image
    print("  Opening image...")
    img = ij.io().open(file_path)
    if img is None:
        print(f"  Failed to open image: {filename}")
        IJ.run("Close All")
        return None

    imp = ij.convert().convert(img, jimport('ij.ImagePlus'))
    if imp is None:
        print(f"  Failed to convert image '{filename}' to ImagePlus.")
        IJ.run("Close All")
        return None

    # Extract fibronectin channel
    print("  Extracting fibronectin channel...")
    if fibronectin_channel > imp.getNChannels() or fibronectin_channel < 1:
        print(f"  Invalid fibronectin channel for image.")
        imp.close()
        IJ.run("Close All")
        return None

    imp_fibronectin = Duplicator().run(
        imp,
        fibronectin_channel,
        fibronectin_channel,
        1,
        imp.getNSlices(),
        1,
        imp.getNFrames()
    )
    imp.close()
    if imp_fibronectin is None:
        print(f"  Failed to extract fibronectin channel.")
        IJ.run("Close All")
        return None
    imp_fibronectin.setTitle(f"{filename}_C{fibronectin_channel}")

    # Reslice
    print("  Performing Reslice...")
    IJ.run(imp_fibronectin, "Reslice [/]...", "output=0.500 start=Top flip rotate avoid")
    imp_fibronectin.close()
    time.sleep(2)

    resliced_imp = IJ.getImage()
    if resliced_imp is None:
        print(f"  Failed to perform Reslice.")
        IJ.run("Close All")
        return None
    resliced_imp.setTitle(f"Reslice_of_{filename}")

    # Z Project
    print("  Performing Z projection...")
    IJ.run(resliced_imp, "Z Project...", "projection=[Max Intensity]")
    resliced_imp.close()
    time.sleep(2)

    projected_imp = IJ.getImage()
    if projected_imp is None:
        print(f"  Failed to perform Z projection.")
        IJ.run("Close All")
        return None
    projected_imp.setTitle(f"MAX_Reslice_of_{filename}")

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

    # Save mask image
    mask_image_path = os.path.join(results_folder, f"Mask_{filename}.tif")
    IJ.saveAs(projected_imp, "Tiff", mask_image_path)
    print(f"  Mask saved to '{mask_image_path}'.")

    # Run Local Thickness
    print("  Running Local Thickness...")
    images_before = set(WindowManager.getImageTitles())
    IJ.runMacro('run("Local Thickness (masked, calibrated, silent)");')
    time.sleep(5)

    images_after = set(WindowManager.getImageTitles())
    new_images = images_after - images_before
    if not new_images:
        print("  Local Thickness plugin failed.")
        IJ.run("Close All")
        return None

    new_image_title = new_images.pop()
    local_thickness_imp = WindowManager.getImage(new_image_title)
    if local_thickness_imp is None:
        print("  Could not retrieve Local Thickness image.")
        IJ.run("Close All")
        return None
    local_thickness_imp.setTitle(f"Local_Thickness_of_{filename}")

    # Set measurements as in macro
    IJ.run("Set Measurements...", "area standard min median redirect=None decimal=3")

    # Clear previous Results
    IJ.run("Clear Results")

    # Measure
    print("  Measuring thickness...")
    IJ.run(local_thickness_imp, "Measure", "")
    rt = jimport('ij.measure.ResultsTable').getResultsTable()

    # Extract results
    if rt is None or rt.getCounter() == 0:
        print("  No measurements.")
        area = std_dev = min_thickness = median_thickness = None
    else:
        try:
            row = rt.getCounter() - 1
            area = rt.getValue("Area", row)
            std_dev = rt.getValue("StdDev", row)
            min_thickness = rt.getValue("Min", row)
            max_thickness = rt.getValue("Max", row)
            median_thickness = rt.getValue("Median", row)
            print(f"  Results - Area: {area}, StdDev: {std_dev}, "
                  f"Min: {min_thickness}, Max: {max_thickness}, Median: {median_thickness}")
        except Exception as e:
            print(f"  Error reading measurements: {e}")
            area = std_dev = min_thickness = max_thickness = median_thickness = None

    # Save thickness image
    thickness_path = os.path.join(results_folder, f"Local_Thickness_{filename}.tif")
    IJ.saveAs(local_thickness_imp, "Tiff", thickness_path)
    print(f"  Local Thickness image saved to '{thickness_path}'.")

    IJ.run("Close All")
    print("  Closed all images.\n")

    return {
        'File_Name': filename,
        'Area': area,
        'StdDev': std_dev,
        'Min': min_thickness,
        'Max': max_thickness,
        'Median': median_thickness
    }


def process_single_folder(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder,
    file_extension,
    fibronectin_channel
):
    """
    Process all image files in a single folder.

    Args:
        ij (imagej.ImageJ): The ImageJ instance.
        IJ, WindowManager, Duplicator, ResultsTable: Java class references.
        folder (str): The folder path to process.
        file_extension (str): The file extension to process (.nd2 or .tiff).
        fibronectin_channel (int): The fibronectin channel index.
    """
    image_files = [
    f for f in os.listdir(folder)
    if f.lower().endswith(file_extension) and not f.startswith(".")
    ]

    if not image_files:
        print(
            f"No '{file_extension}' files found in folder '{folder}'. "
            "Skipping."
        )
        return

    print(f"\nProcessing folder: {folder}")
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    results_folder_name = f"Thickness_assay_results_{timestamp}"
    results_folder = os.path.join(folder, results_folder_name)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    print(f"Results will be saved in: {results_folder}")

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
            file_extension,
            results_folder
        )
        if result is not None:
            summary_data.append(result)

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file_path = os.path.join(results_folder, 'Thickness_Summary.csv')
        summary_df.to_csv(summary_file_path, index=False)
        print(f"Folder analysis complete. Data saved to '{summary_file_path}'.")
    else:
        print(f"No data to save in summary for folder '{folder}'.")


def process_all_folders(
    ij,
    IJ,
    WindowManager,
    Duplicator,
    ResultsTable,
    folder_paths,
    file_extension,
    fibronectin_channel
):
    """
    Process all provided folders.

    Args:
        ij (imagej.ImageJ): The ImageJ instance.
        IJ, WindowManager, Duplicator, ResultsTable: Java class references.
        folder_paths (list[str]): List of folder paths to process.
        file_extension (str): The file extension to process (.nd2 or .tiff).
        fibronectin_channel (int): The fibronectin channel index.
    """
    for folder in folder_paths:
        process_single_folder(
            ij,
            IJ,
            WindowManager,
            Duplicator,
            ResultsTable,
            folder,
            file_extension,
            fibronectin_channel
        )
    print("\nAll folders have been processed.")


def main():
    """
    Main execution function.
    """
    parser = argparse.ArgumentParser(description='Thickness analysis script')
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help='Path to the input JSON file containing folder paths'
    )
    args = parser.parse_args()

    print("Initializing ImageJ...")
    ij = initialize_imagej()
    print("ImageJ initialization completed.")

    (
        IJ,
        ImagePlus,
        WindowManager,
        ResultsTable,
        Duplicator,
        System
    ) = import_java_classes()

    folder_paths = get_folder_paths(args.input)
    file_extension = get_file_type_choice()
    num_channels = get_num_channels()
    fibronectin_channel = get_fibronectin_channel(num_channels)

    if not confirm_processing():
        print("Processing aborted by user.")
        ij.dispose()
        System.exit(0)

    process_all_folders(
        ij,
        IJ,
        WindowManager,
        Duplicator,
        ResultsTable,
        folder_paths,
        file_extension,
        fibronectin_channel
    )

    print("Disposing of ImageJ context...")
    ij.dispose()
    print("Exiting the program.")
    System.exit(0)


if __name__ == "__main__":
    main()
