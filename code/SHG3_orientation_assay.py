#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Orientation Analysis Script for Fibers (Single-Line Summary Per Image)

This script processes images (e.g., .tif, .nd2) to analyze fiber orientation and measure total area
via the ParticleAnalyzer summary, generating exactly one row per image in the final CSV.

Main steps:
1. Creates 2D projections (max-intensity) from 3D images.
2. Thresholds a duplicate and runs ParticleAnalyzer in summary mode => one summary row per image.
3. Performs gradient-based orientation analysis (using orientationpy) on the original 2D projection.
4. Summarizes the distribution of angles and saves CSV/plots.

All comments/documentation are in English.
"""

import imagej
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
import sys
import numpy as np
from skimage import io
import scyjava as sj
import argparse
import json
import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# orientationpy library for gradient-based orientation
from orientationpy import (
    computeGradient,
    computeStructureTensor,
    computeIntensity,
    computeStructureDirectionality,
    computeOrientation,
)

# Set a non-interactive backend for matplotlib
matplotlib.use('Agg')


class ImageJInitializationError(Exception):
    """Exception raised for unsuccessful initialization of ImageJ."""
    pass

class LifFileProcessorError(Exception):
    """Custom exception for LIF file processing errors (if applicable)."""
    pass


def initialize_imagej():
    """
    Initializes the ImageJ context in headless mode with Bio-Formats support.
    Returns the initialized ImageJ context.
    Raises ImageJInitializationError if ImageJ fails to initialize.
    """
    print("Initializing ImageJ in headless mode...")
    try:
        # Use mode='headless'
        ij_instance = imagej.init('sc.fiji:fiji', mode='headless')
    except Exception as e:
        raise ImageJInitializationError(f"Failed to initialize ImageJ: {e}")
    print("ImageJ successfully initialized.")
    return ij_instance


def get_paths_to_files(input_file_path):
    """
    Reads an input JSON file containing folder paths under the key 'paths_to_files'.

    Returns a list of valid folder paths. Prints info about each folder's contents.
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"File '{input_file_path}' does not exist.")

    if not input_file_path.lower().endswith('.json'):
        raise ValueError("Input file must be in .json format.")

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
            print(f"\nFolder: {folder_path}\n  Number of files: {num_files}\n  File types: {', '.join(exts)}")
            valid_paths_to_files.append(folder_path)
        else:
            print(f"\nFolder '{folder_path}' does not exist.")

    if not valid_paths_to_files:
        raise ValueError("No valid folders found in 'paths_to_files'.")

    print(f"\nFound {len(valid_paths_to_files)} folders for processing.")
    return valid_paths_to_files


def get_angle_value():
    """Prompt user for angle. Default=15 if blank or invalid."""
    angle_input = input("\nEnter angle for analysis (default is 15): ").strip()
    if angle_input == '':
        return 15.0
    else:
        try:
            return float(angle_input)
        except ValueError:
            print("Invalid input. Using default=15.")
            return 15.0


def get_fiber_channel_indices(paths_to_files):
    """
    Prompts user for fiber channel index (1-based) for each folder or one global index if user chooses "yes".
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


def create_results_folders(folder_path, angle_value_str, timestamp):
    """
    Create main results folder: 'Alignment_assay_results_angle_{angle_value_str}_{timestamp}'
    plus subfolders 'Tables' and 'Images'.
    Returns (results_folder, table_folder, images_folder).
    """
    results_dirname = f"Alignment_assay_results_angle_{angle_value_str}_{timestamp}"
    results_folder = os.path.join(folder_path, results_dirname)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    print(f"Results folder: {results_folder}")

    table_folder = os.path.join(results_folder, "Tables")
    images_folder = os.path.join(results_folder, "Images")
    Path(table_folder).mkdir(exist_ok=True)
    Path(images_folder).mkdir(exist_ok=True)
    print(f"  Created subfolders: 'Tables', 'Images'")

    return results_folder, table_folder, images_folder


def analyze_particles_summary(imp, mask_imp, min_size=0.0, max_size=1e9, min_circ=0.0, max_circ=1.0):
    """
    Analyze particles on 'imp' using 'mask_imp' in SUMMARY mode (one row per image).
    Measures AREA, MEDIAN INTENSITY, and other parameters.

    Returns:
        sumRT: The summary ResultsTable (one row).
    """
    ParticleAnalyzer = sj.jimport('ij.plugin.filter.ParticleAnalyzer')
    ResultsTable = sj.jimport('ij.measure.ResultsTable')
    Measurements = sj.jimport('ij.measure.Measurements')
    IJ = sj.jimport('ij.IJ')
    ImageProcessor = sj.jimport('ij.process.ImageProcessor')

    # Combine flags for AREA, MEDIAN, etc.
    measurements = (Measurements.AREA |
                   Measurements.MEDIAN |
                   Measurements.MEAN |
                   Measurements.INTEGRATED_DENSITY |
                   Measurements.LIMIT)
    # Show summary only (no per-particle table)
    options = ParticleAnalyzer.SHOW_SUMMARY

    # We'll create a separate table for the summary
    summaryRT = ResultsTable()
    # By default, ParticleAnalyzer uses a static Summary table. We override that:
    ParticleAnalyzer.setSummaryTable(summaryRT)

    # This table is for per-particle data (we won't use it)
    rt = ResultsTable()

    # Get the image processors
    imp_processor = imp.getProcessor()
    mask_processor = mask_imp.getProcessor()

    # Create a new image processor with the mask applied
    masked_processor = imp_processor.duplicate()
    masked_processor.copyBits(mask_processor, 0, 0, ImageProcessor.BLACK)

    # Create a new image with the masked processor
    masked_imp = sj.jimport('ij.ImagePlus')("Masked Image", masked_processor)

    # Ensure the image is thresholded
    IJ.setAutoThreshold(masked_imp, "Moments")
    Prefs = sj.jimport('ij.Prefs')
    Prefs.blackBackground = True
    IJ.run(masked_imp, "Convert to Mask", "")

    # Analyze the masked image
    pa = ParticleAnalyzer(options, measurements, rt,
                          min_size, max_size,
                          min_circ, max_circ)

    success = pa.analyze(masked_imp)
    if not success:
        raise ValueError("ParticleAnalysis failed (summary).")

    # The summary table now has a single row for this image
    return summaryRT

def process_part1(folder_path, results_folder, fiber_channel_index, desired_width, desired_height, ij):
    """
    Part 1: 
    1) For each image, create a 2D projection (max-intensity).
    2) Save it as <filename>_processed.tif (un-thresholded).
    3) Duplicate that projection => threshold => show the mask => run 'analyze_particles_summary'
    4) We gather only the single-line summary for each file into a single CSV.

    Returns a dict mapping processed_base_name -> {some info} for part3.
    """
    print("\n[Part 1] Creating 2D projections and single-line area summaries...")

    IJ = sj.jimport('ij.IJ')
    ZProjector = sj.jimport('ij.plugin.ZProjector')
    Prefs = sj.jimport('ij.Prefs')

    # Folder for storing threshold masks
    mask_folder = os.path.join(results_folder, "Masks")
    os.makedirs(mask_folder, exist_ok=True)

    all_summaries = []
    z_stacks_info_folder = {}

    IJ.run("Close All")

    for filename in os.listdir(folder_path):
        # Only .tif/.tiff/.nd2
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
            print(f"  Fiber channel {fiber_channel_index} out of range for '{filename}'. Skipping.")
            imp.close()
            continue

        # Duplicate the fiber channel
        imp.setC(fiber_channel_index)
        IJ.run(imp, "Duplicate...", f"title=imp_fiber duplicate channels={fiber_channel_index}")
        imp_fiber = IJ.getImage()

        # Max-intensity projection
        zp = ZProjector(imp_fiber)
        zp.setMethod(ZProjector.MAX_METHOD)
        zp.doProjection()
        fiber_proj = zp.getProjection()

        # Resize and convert to 16-bit
        fiber_proj = fiber_proj.resize(desired_width, desired_height, "bilinear")
        # IJ.run(fiber_proj, "16-bit", "")

        # Enhance contrast
        # IJ.run(fiber_proj, "Enhance Contrast...", "saturated=0.35")

        # Save the un-thresholded projection 
        out_basename = os.path.splitext(filename)[0] + "_processed.tif"
        out_path = os.path.join(results_folder, out_basename)
        IJ.saveAs(fiber_proj, "Tiff", out_path)
        print(f"  Saved 2D projection => {out_path}")

        # Create a copy for thresholding
        IJ.run("Close All")
        area_imp = IJ.openImage(out_path)
        if area_imp is None:
            print(f"  Could not reopen '{out_path}' for thresholding.")
            fiber_proj.close()
            imp_fiber.close()
            imp.close()
            continue

        # Threshold => mask
        IJ.setAutoThreshold(area_imp, "Moments")
        Prefs.blackBackground = True
        IJ.run(area_imp, "Convert to Mask", "")

        # Save the thresholded mask
        mask_name = os.path.splitext(filename)[0] + "_mask.tif"
        mask_path = os.path.join(mask_folder, mask_name)
        IJ.saveAs(area_imp, "Tiff", mask_path)
        print(f"  Saved mask => {mask_path}")

        # Analyze the original image using the mask
        try:
            sumRT = analyze_particles_summary(fiber_proj, area_imp)
        except ValueError as e:
            print(f"  Particle analysis failed => {e}")
            area_imp.close()
            fiber_proj.close()
            imp_fiber.close()
            imp.close()
            continue

        # The summary table has 1 row with headings like: 
        # ['Slice', 'Count', 'Total Area', 'Average Size', '%Area', 'Avg. Circ.', 'Mean Gray Value']
        headings = sumRT.getHeadings()  
        # Usually only 1 row in the summary
        if sumRT.getCounter() < 1:
            print("  Summary table is empty?!")
        else:
            row_idx = 0
            col_data = {}
            for col_name in headings:
                val = sumRT.getValue(col_name, row_idx)
                col_data[col_name] = [val]

            # Make a DataFrame with one row
            df_summary = pd.DataFrame(col_data)
            df_summary["Image_Filename"] = [filename]
            all_summaries.append(df_summary)

        # Cleanup
        area_imp.close()
        fiber_proj.close()
        imp_fiber.close()
        imp.close()
        IJ.run("Close All")

        # For part3
        processed_name_no_ext = os.path.splitext(out_basename)[0]
        z_stacks_info_folder[processed_name_no_ext] = {
            "original_filename": filename,
            "number_of_z_stacks": slices,
            "z_stack_type": "slices"
        }

    # Combine all summaries => single CSV
    if all_summaries:
        df_all = pd.concat(all_summaries, ignore_index=True)
        all_area_csv = os.path.join(results_folder, "All_area_measurements.csv")
        df_all.to_csv(all_area_csv, index=False)
        print(f"\nSaved single-line summaries to => {all_area_csv}")
    else:
        print("\nNo summary data found. Possibly no images processed?")

    print(f"\nPart 1 complete for '{folder_path}'.")
    return z_stacks_info_folder

def process_part2_orientation_analysis(results_folder, images_folder):
    """
    Part 2: Perform orientation analysis on the un-thresholded 2D-projected images (*_processed.tif).
    Saves orientation distributions and overlay images.
    """
    print("\n[Part 2] Orientation Analysis...")

    processed_files = [f for f in os.listdir(results_folder) if f.lower().endswith('_processed.tif')]
    if not processed_files:
        print("No processed TIFs found for orientation analysis.")
        return

    for fname in processed_files:
        fp = os.path.join(results_folder, fname)
        print(f"\nAnalyzing orientation => {fp}")

        image_gray = io.imread(fp).astype(float)
        if image_gray.size == 0:
            print("Empty image? Skipping.")
            continue

        # Compute gradients
        anisotropy = np.array([1., 1., 1.])
        Gy, Gx = computeGradient(image_gray, mode='splines', anisotropy=anisotropy)
        # Structure tensor
        st = computeStructureTensor((Gy, Gx), sigma=1)
        directionality = computeStructureDirectionality(st)
        orientations = computeOrientation(st)

        # Build orientation histogram
        orient_flat = orientations["theta"].flatten()
        orient_flat = orient_flat[~np.isnan(orient_flat)]
        hist, bin_edges = np.histogram(orient_flat, bins=180, range=(-90, 90))
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Save CSV
        tables_dir = os.path.join(results_folder, "Tables")
        os.makedirs(tables_dir, exist_ok=True)
        out_csv = os.path.splitext(fname)[0] + "_orientation_distribution.csv"
        out_csv_fp = os.path.join(tables_dir, out_csv)
        df = pd.DataFrame({"orientation_angle": bin_centers, "occurrence_value": hist})
        df.to_csv(out_csv_fp, index=False)
        print(f"  Orientation CSV => {out_csv_fp}")

        # HSV overlay
        vmin, vmax = 10, 1e8
        clipped_dir = np.clip(directionality, vmin, vmax)
        clipped_log = np.log(clipped_dir)
        clipped_log -= clipped_log.min()
        clipped_log /= clipped_log.max()
        clipped_log[image_gray == 0] = 0

        imHsv = np.zeros((image_gray.shape[0], image_gray.shape[1], 3), dtype="float32")
        imHsv[:, :, 0] = (orientations["theta"] + 90) / 180.0
        imHsv[:, :, 1] = clipped_log
        imHsv[:, :, 2] = image_gray / (image_gray.max() if image_gray.max() > 0 else 1)

        fig, ax = plt.subplots(figsize=(6,6))
        rgb = matplotlib.colors.hsv_to_rgb(imHsv)
        ax.imshow(rgb)
        ax.set_title(f"Orientation Overlay: {fname}")

        sm = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-90, vmax=90), cmap='hsv')
        sm.set_array([])
        fig.colorbar(sm, ax=ax, orientation="vertical", label="Angle (degrees)")

        overlay_name = os.path.splitext(fname)[0] + "_orientation_overlay.png"
        overlay_path = os.path.join(images_folder, overlay_name)
        fig.savefig(overlay_path)
        plt.close(fig)
        print(f"  Orientation overlay => {overlay_path}")


def process_part3(results_folder, analysis_folder, angle_value, z_stacks_info):
    """
    Part 3: Summarize orientation distribution results.

    Reads each *_orientation_distribution.csv, re-centers angles around the max-occurrence angle,
    determines the percentage of angles within +/- angle_value => alignment vs disorganized.
    Saves a final summary "Alignment_Summary.csv".
    """
    print("\n[Part 3] Summarizing orientation distributions...")

    tables_dir = os.path.join(results_folder, "Tables")
    if not os.path.isdir(tables_dir):
        print("No 'Tables' folder found, orientation analysis might be skipped.")
        return

    csv_files = [f for f in os.listdir(tables_dir) if f.endswith('_orientation_distribution.csv')]
    if not csv_files:
        print("No orientation CSV files found.")
        return

    summary_data = []

    for fname in csv_files:
        fp = os.path.join(tables_dir, fname)
        print(f"\nProcessing distribution => {fp}")

        df_in = pd.read_csv(fp)
        df_in.columns = ["orientation_angle", "occurrence_value"]

        max_occ_val = df_in["occurrence_value"].max()
        angle_of_max_occ = df_in["orientation_angle"][df_in["occurrence_value"].idxmax()]

        df_in["angles_normalized_to_angle_of_MOV"] = df_in["orientation_angle"] - angle_of_max_occ
        def correct_angle(x):
            if x < -90:
                return x + 180
            elif x > 90:
                return x - 180
            else:
                return x
        df_in["corrected_angles"] = df_in["angles_normalized_to_angle_of_MOV"].apply(correct_angle)
        df_in["rank_of_angle_occurrence_value"] = df_in["corrected_angles"].rank(method='min')

        total_occ = df_in["occurrence_value"].sum()
        df_in["percent_occurrence_value_to_sum_of_occurrence_value"] = (df_in["occurrence_value"] / total_occ) * 100

        # filter +/- angle_value
        f_data = df_in[
            (df_in["corrected_angles"] >= -angle_value) &
            (df_in["corrected_angles"] <= angle_value)
        ]
        aligned_perc = f_data["percent_occurrence_value_to_sum_of_occurrence_value"].sum()
        orientation_mode = "aligned" if aligned_perc >= 55 else "disorganized"

        df_sorted = df_in.sort_values(by="rank_of_angle_occurrence_value")

        base_name = os.path.splitext(fname.replace('_orientation_distribution.csv', ''))[0]
        out_proc_name = base_name + "_processed.csv"
        out_proc_path = os.path.join(analysis_folder, out_proc_name)
        df_sorted.to_csv(out_proc_path, index=False)
        print(f"  Processed orientation => {out_proc_path}")

        # Info from part1
        z_info = z_stacks_info.get(base_name, {})
        stacks_count = z_info.get("number_of_z_stacks", "N/A")
        stack_type = z_info.get("z_stack_type", "N/A")

        summary_data.append({
            "File_Name": fname,
            "Number_of_Z_Stacks": stacks_count,
            "Z_Stack_Type": stack_type,
            f"Percentage_Fibers_Aligned_Within_{angle_value}_Degree": aligned_perc,
            "Orientation_Mode": orientation_mode
        })

    df_summary = pd.DataFrame(summary_data)
    align_sum = os.path.join(analysis_folder, "Alignment_Summary.csv")
    df_summary.to_csv(align_sum, index=False)
    print(f"\nSummary => {align_sum}\nPart 3 done.")


def process_folder(folder_path, fiber_channel_index, angle_value, desired_width, desired_height, ij):
    """
    Processes a single folder in 3 steps:
      1) 2D projection + threshold => one-line summary from ParticleAnalyzer
      2) Orientation analysis on un-thresholded image
      3) Summarize orientation CSV => final alignment summary
    """
    if not os.path.isdir(folder_path):
        print(f"Folder '{folder_path}' does not exist, skipping.")
        return

    angle_str = str(angle_value).replace('.', '_')
    now = datetime.now().strftime("%Y%m%d_%H%M%S")

    results_folder, table_folder, images_folder = create_results_folders(folder_path, angle_str, now)

    # Part 1: Summaries
    z_stacks_info = process_part1(folder_path, results_folder,
                                  fiber_channel_index,
                                  desired_width, desired_height, ij)

    # Part 2: Orientation
    process_part2_orientation_analysis(results_folder, images_folder)

    # Part 3: Summarize orientation
    analysis_folder = os.path.join(results_folder, "Analysis")
    Path(analysis_folder).mkdir(exist_ok=True)
    process_part3(results_folder, analysis_folder, angle_value, z_stacks_info)


def main():
    parser = argparse.ArgumentParser(description="Fiber orientation & single-line summary area script")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file with 'paths_to_files'")
    args = parser.parse_args()

    # 1) Initialize ImageJ
    try:
        ij = initialize_imagej()
    except ImageJInitializationError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # 2) Desired output size
    desired_width = 1024
    desired_height = 1024

    # 3) Paths to folders from JSON
    try:
        folder_paths = get_paths_to_files(args.input)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading input paths: {e}")
        sys.exit(1)

    # 4) Angle
    angle_val = get_angle_value()

    # 5) Fiber channel indices
    fiber_channels = get_fiber_channel_indices(folder_paths)

    # Confirm
    ans = input("\nStart processing? (yes/no): ").lower().strip()
    if ans not in ("yes", "y"):
        print("Canceled by user.")
        sys.exit(0)

    # 6) Process each folder
    for fpath in folder_paths:
        fib_ch = fiber_channels.get(fpath, 1)
        process_folder(fpath, fib_ch, angle_val, desired_width, desired_height, ij)

    print("\nAll done.")

    # 7) Dispose ImageJ
    print("Closing ImageJ context...")
    ij.context().dispose()
    print("Script finished.")


if __name__ == "__main__":
    main()