#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mask Creation & Intensity/Area Measurement Script (Single-Line Summary Per Image)

Key Features:
1. User can select auto threshold or manual threshold (enter lower and upper).
2. For each folder listed in a JSON, user can specify which fiber channel to use.
3. If an image is a Z-stack, the script creates a max-intensity projection on the chosen channel.
4. It saves the projected image as the "original" for further measurements.
5. Then it creates a binary mask (0/255) by thresholding that projection.
6. A bitwise AND operation merges the original intensities with the mask, leaving only
   pixels within the mask as non-zero.
7. The script saves both the mask (binary) and the masked image (intensities with background=0).
8. It performs `ParticleAnalyzer` in summary mode to measure Area, Mean, Median, IntDen, etc.
9. It writes a CSV (All_area_measurements.csv) summarizing these measurements for each image.
"""

import imagej
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
import sys
import argparse
import json
import scyjava as sj


class ImageJInitializationError(Exception):
    """Exception raised for unsuccessful initialization of ImageJ."""
    pass


def initialize_imagej():
    """
    Initializes the ImageJ context in headless mode with Bio-Formats support.
    Returns the initialized ImageJ context.
    Raises ImageJInitializationError if ImageJ fails to initialize.
    """
    print("Initializing ImageJ in headless mode...")
    try:
        ij_instance = imagej.init('sc.fiji:fiji', mode='headless')
    except Exception as e:
        raise ImageJInitializationError(f"Failed to initialize ImageJ: {e}")
    print("ImageJ successfully initialized.")
    return ij_instance


def get_paths_to_files(input_file_path):
    """
    Reads an input JSON file containing folder paths under the key 'paths_to_files'.
    Returns a list of valid folder paths. Prints basic info about each folder's contents.
    
    The JSON format should be something like:
    
    {
      "paths_to_files": [
        "/path/to/folderA",
        "/path/to/folderB"
      ]
    }
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"File '{input_file_path}' does not exist.")

    if not input_file_path.lower().endswith('.json'):
        raise ValueError("Input file must be a .json file.")

    with open(input_file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
        paths_to_files = data.get('paths_to_files', [])

    if not paths_to_files:
        raise ValueError("JSON does not contain 'paths_to_files' or it is empty.")

    valid_paths_to_files = []
    for folder_path in paths_to_files:
        if os.path.isdir(folder_path):
            files = os.listdir(folder_path)
            num_files = len(files)
            exts = set([os.path.splitext(f)[1].lower() for f in files])
            print(f"\nFolder: {folder_path}")
            print(f"  Number of files: {num_files}")
            print(f"  File types: {', '.join(exts)}")
            valid_paths_to_files.append(folder_path)
        else:
            print(f"\nFolder '{folder_path}' does not exist or is invalid.")

    if not valid_paths_to_files:
        raise ValueError("No valid folders found in 'paths_to_files'.")

    print(f"\nFound {len(valid_paths_to_files)} valid folder(s) for processing.")
    return valid_paths_to_files


def get_fiber_channel_indices(paths_to_files):
    """
    Prompts user for the fiber channel index (1-based) for each folder
    or one global index if user chooses "yes".
    Returns a dict mapping folder_path -> channel_index.
    """
    fiber_channel_indices = {}
    same_channel = input("\nDo all folders use the same fiber channel? (yes/no): ").strip().lower()
    if same_channel in ('yes', 'y'):
        ch_str = input("Enter fiber channel index (1-based): ").strip()
        if ch_str.isdigit():
            fib_ch = int(ch_str)
        else:
            print("Invalid input. Using default channel=1.")
            fib_ch = 1
        for folder in paths_to_files:
            fiber_channel_indices[folder] = fib_ch
    else:
        for folder in paths_to_files:
            print(f"\nFolder: {folder}")
            ch_str = input("Enter fiber channel index (1-based): ").strip()
            if ch_str.isdigit():
                fib_ch = int(ch_str)
            else:
                print("Invalid input. Using default channel=1.")
                fib_ch = 1
            fiber_channel_indices[folder] = fib_ch

    return fiber_channel_indices


def get_threshold_choice():
    """
    Asks the user whether to use 'auto' threshold or 'manual' threshold.
    If 'manual', user enters both lower and upper threshold.
    Returns a tuple: (mode, (lower, upper)) or (mode, None).
    """
    mode = input("\nUse auto threshold or manual threshold? (auto/manual): ").strip().lower()
    if mode == 'manual':
        # Ask for lower threshold
        lower_str = input("Enter lower threshold (e.g. 1000): ").strip()
        try:
            lower_val = float(lower_str)
        except ValueError:
            print("Invalid input. Using default lower=0.")
            lower_val = 0

        # Ask for upper threshold
        upper_str = input("Enter upper threshold (e.g. 65535): ").strip()
        try:
            upper_val = float(upper_str)
        except ValueError:
            print("Invalid input. Using default upper=65535.")
            upper_val = 65535

        return "manual", (lower_val, upper_val)
    else:
        # Default to auto
        return "auto", None


def create_results_folders(folder_path, timestamp):
    """
    Creates a main results folder: 'Mask_Intensity_Area_Results_{timestamp}'
    plus subfolders: 'Tables', 'Masks', and 'Masked_Images'.

    Returns (results_folder, table_folder, mask_folder, masked_folder).
    """
    results_dirname = f"Mask_Intensity_Area_Results_{timestamp}"
    results_folder = os.path.join(folder_path, results_dirname)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    print(f"Results folder: {results_folder}")

    table_folder = os.path.join(results_folder, "Tables")
    Path(table_folder).mkdir(exist_ok=True)

    mask_folder = os.path.join(results_folder, "Masks")
    Path(mask_folder).mkdir(exist_ok=True)

    masked_folder = os.path.join(results_folder, "Masked_Images")
    Path(masked_folder).mkdir(exist_ok=True)

    print("  Subfolders created:")
    print(f"    -> {table_folder}")
    print(f"    -> {mask_folder}")
    print(f"    -> {masked_folder}")

    return results_folder, table_folder, mask_folder, masked_folder


def analyze_particles_summary(imp_original, imp_mask, masked_out_path,
                              min_size=0.0, max_size=1e9,
                              min_circ=0.0, max_circ=1.0):
    """
    Analyze particles on the 'original' image by bitwise AND with a (binary) mask,
    then measure area, mean, median, integrated density, etc.
    Saves a smoothed + contrast-enhanced copy of the masked image but uses the
    original masked image for analysis.
    """
    # Java imports
    ParticleAnalyzer = sj.jimport('ij.plugin.filter.ParticleAnalyzer')
    ResultsTable = sj.jimport('ij.measure.ResultsTable')
    Measurements = sj.jimport('ij.measure.Measurements')
    IJ = sj.jimport('ij.IJ')
    Prefs = sj.jimport('ij.Prefs')
    Blitter = sj.jimport('ij.process.Blitter')

    # Combine flags for AREA, MEDIAN, MEAN, INT_DEN, etc.
    measurements = (Measurements.AREA |
                    Measurements.MEDIAN |
                    Measurements.MEAN |
                    Measurements.INTEGRATED_DENSITY |
                    Measurements.LIMIT)
    # Show summary only
    options = ParticleAnalyzer.SHOW_SUMMARY

    # Create a summary table
    summaryRT = ResultsTable()
    ParticleAnalyzer.setSummaryTable(summaryRT)

    # 1) Duplicate the original image processor
    imp_processor = imp_original.getProcessor().duplicate()
    mask_processor = imp_mask.getProcessor()

    # 2) Bitwise AND -> background=0, object intensities retained
    imp_processor.copyBits(mask_processor, 0, 0, Blitter.AND)

    # This masked_imp is for actual measurements
    masked_imp = sj.jimport('ij.ImagePlus')("Masked Image", imp_processor)

    # --- For visualization only: smooth + enhance contrast + save ---
    masked_imp_for_save = masked_imp.duplicate()
    # Apply smoothing twice
    IJ.run(masked_imp_for_save, "Smooth", "")
    IJ.run(masked_imp_for_save, "Smooth", "")
    # Enhance contrast
    IJ.run(masked_imp_for_save, "Enhance Contrast...", "saturated=0.35")
    # Save the processed image (visualization copy)
    IJ.saveAs(masked_imp_for_save, "Tiff", masked_out_path)
    # Close the visualization copy
    masked_imp_for_save.close()

    # --- Now threshold the analysis version from 1..∞ for ParticleAnalyzer
    IJ.setThreshold(masked_imp, 1, float('inf'))
    Prefs.blackBackground = True

    # ParticleAnalyzer in summary mode
    pa = ParticleAnalyzer(options, measurements, ResultsTable(),
                          min_size, max_size,
                          min_circ, max_circ)
    success = pa.analyze(masked_imp)
    if not success:
        raise ValueError("ParticleAnalyzer failed in summary mode.")

    return summaryRT


def process_images_in_folder(folder_path, fiber_channel_index, ij,
                             threshold_mode, manual_threshold_range,
                             desired_width=1024, desired_height=1024):
    """
    Core logic for one folder:
      1) For each image, if it has multiple Z slices, do a Max Intensity Projection on the fiber channel.
      2) Save that 2D projection as <filename>_processed.tif.
      3) Make a copy -> threshold (auto or manual) -> Convert to Mask (0/255).
      4) Save the binary mask.
      5) Merge the original with the mask by bitwise AND => saves the masked image.
      6) Run ParticleAnalyzer in summary mode -> get one CSV row per image.
      7) Collect all results into "All_area_measurements.csv".
    """
    if not os.path.isdir(folder_path):
        print(f"Folder '{folder_path}' does not exist, skipping.")
        return

    # Java classes
    IJ = sj.jimport('ij.IJ')
    ZProjector = sj.jimport('ij.plugin.ZProjector')
    Prefs = sj.jimport('ij.Prefs')

    # Create the main results folders
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_folder, table_folder, mask_folder, masked_folder = create_results_folders(folder_path, now)

    all_summaries = []

    # Close all images to avoid confusion
    IJ.run("Close All")

    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(('.tif', '.tiff', '.nd2')):
            print(f"Skipping '{filename}' (not .tif/.tiff/.nd2).")
            continue

        file_path = os.path.join(folder_path, filename)
        print(f"\nProcessing: {file_path}")

        IJ.run("Close All")
        imp = IJ.openImage(file_path)
        if imp is None:
            print(f"Failed to open '{file_path}'.")
            continue

        w, h, channels, slices, frames = imp.getDimensions()
        print(f"  Dimensions => w={w}, h={h}, channels={channels}, slices={slices}, frames={frames}")

        if fiber_channel_index > channels:
            print(f"  Fiber channel {fiber_channel_index} is out of range for '{filename}'. Skipping.")
            imp.close()
            continue

        # Select the desired fiber channel
        imp.setC(fiber_channel_index)
        # Duplicate only this channel
        IJ.run(imp, "Duplicate...", f"title=imp_fiber duplicate channels={fiber_channel_index}")
        imp_fiber = IJ.getImage()

        # If multiple Z slices, do a Max Intensity Projection
        if slices > 1:
            zp = ZProjector(imp_fiber)
            zp.setMethod(ZProjector.MAX_METHOD)
            zp.doProjection()
            fiber_proj = zp.getProjection()
        else:
            fiber_proj = imp_fiber.duplicate()

        # Resize if desired
        fiber_proj = fiber_proj.resize(desired_width, desired_height, "bilinear")

        # Save this 2D projection as the "original" for intensity analysis
        out_basename = os.path.splitext(filename)[0] + "_processed.tif"
        out_path = os.path.join(results_folder, out_basename)
        IJ.saveAs(fiber_proj, "Tiff", out_path)
        print(f"  Saved 2D projection => {out_path}")

        # Create the binary mask
        IJ.run("Close All")
        imp_for_mask = IJ.openImage(out_path)
        if imp_for_mask is None:
            print(f"  Could not reopen '{out_path}' for thresholding.")
            fiber_proj.close()
            imp_fiber.close()
            imp.close()
            continue

        Prefs.blackBackground = True

        # Depending on user choice, set auto or manual threshold
        if threshold_mode == "manual":
            lower_val, upper_val = manual_threshold_range
            print(f"  Using manual threshold => [{lower_val}, {upper_val}]")
            IJ.setThreshold(imp_for_mask, lower_val, upper_val)
            IJ.run(imp_for_mask, "Convert to Mask", "")
        else:
            print("  Using auto threshold (Moments).")
            IJ.setAutoThreshold(imp_for_mask, "Moments")
            IJ.run(imp_for_mask, "Convert to Mask", "")

        # Optionally invert the mask if needed:
        # IJ.run(imp_for_mask, "Invert", "")

        # Save the 0/255 binary mask
        mask_name = os.path.splitext(filename)[0] + "_mask.tif"
        mask_path = os.path.join(mask_folder, mask_name)
        IJ.saveAs(imp_for_mask, "Tiff", mask_path)
        print(f"  Saved binary mask => {mask_path}")

        # Reopen the "original" image to measure from
        IJ.run("Close All")
        imp_original = IJ.openImage(out_path)
        if imp_original is None:
            print(f"  Could not reopen '{out_path}' as imp_original.")
            imp_for_mask.close()
            fiber_proj.close()
            imp_fiber.close()
            imp.close()
            continue

        # We'll also specify the final masked image filename
        masked_image_name = os.path.splitext(filename)[0] + "_masked.tif"
        masked_out_path = os.path.join(masked_folder, masked_image_name)

        # Perform Particle Analysis with the mask
        try:
            sumRT = analyze_particles_summary(imp_original, imp_for_mask, masked_out_path)
        except ValueError as e:
            print(f"  Particle analysis failed => {e}")
            imp_for_mask.close()
            fiber_proj.close()
            imp_fiber.close()
            imp.close()
            continue

        # Extract the single summary row
        headings = sumRT.getHeadings()
        if sumRT.getCounter() < 1:
            print("  Summary table is empty?!")
        else:
            row_idx = 0
            col_data = {}
            for col_name in headings:
                val = sumRT.getValue(col_name, row_idx)
                col_data[col_name] = [val]

            # Convert to DataFrame for convenience
            df_summary = pd.DataFrame(col_data)
            df_summary["Image_Filename"] = [filename]
            all_summaries.append(df_summary)

        # Cleanup
        imp_for_mask.close()
        fiber_proj.close()
        imp_fiber.close()
        imp_original.close()
        imp.close()
        IJ.run("Close All")

    # After processing all images in this folder, combine summary rows
    if all_summaries:
        df_all = pd.concat(all_summaries, ignore_index=True)
        all_area_csv = os.path.join(table_folder, "All_area_measurements.csv")
        df_all.to_csv(all_area_csv, index=False)
        print(f"\nSingle-line summaries saved to => {all_area_csv}")
    else:
        print("\nNo summary data found. Possibly no images processed?")

    print(f"\nFinished processing folder: '{folder_path}'.\n")


def main():
    parser = argparse.ArgumentParser(description="Mask + AND with original -> measure Intensity/Area (auto/manual threshold).")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file with 'paths_to_files'")
    args = parser.parse_args()

    # 1) Initialize ImageJ
    try:
        ij = initialize_imagej()
    except ImageJInitializationError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # 2) Read folders from JSON
    try:
        folder_paths = get_paths_to_files(args.input)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading input paths: {e}")
        sys.exit(1)

    # 3) Ask for fiber channel(s)
    fiber_channels = get_fiber_channel_indices(folder_paths)

    # 4) Ask user if threshold is auto or manual; if manual, get lower/upper
    threshold_mode, manual_range = get_threshold_choice()

    # 5) Confirm
    ans = input("\nStart processing? (yes/no): ").lower().strip()
    if ans not in ("yes", "y"):
        print("Canceled by user.")
        sys.exit(0)

    # 6) Process each folder
    for fpath in folder_paths:
        fib_ch = fiber_channels.get(fpath, 1)
        process_images_in_folder(
            folder_path=fpath,
            fiber_channel_index=fib_ch,
            ij=ij,
            threshold_mode=threshold_mode,
            manual_threshold_range=manual_range
        )

    print("All done.")

    # 7) Dispose the ImageJ context
    print("Closing ImageJ context...")
    ij.context().dispose()
    print("Script finished.")


if __name__ == "__main__":
    main()
