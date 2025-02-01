#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import json
import os
from datetime import datetime
from pathlib import Path

import imagej
import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scyjava as sj
from skimage import io

from orientationpy import (computeGradient,
                           computeOrientation,
                           computeStructureDirectionality,
                           computeStructureTensor)


# Custom exception classes
class ImageJInitializationError(Exception):
    """Exception raised for unsuccessful
    initialization of ImageJ."""
    pass


# Helper function: apply angle correction
def correct_angle(x):
    if x < -90:
        return x + 180
    elif x > 90:
        return x - 180
    else:
        return x


def get_folder_paths(input_file_path):
    """
    Reads an input JSON file containing folder paths.

    Args:
        input_file_path (str): Path to the input JSON file.

    Returns:
        List[str]: List of folder paths.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file does not contain folder paths.
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(f"File '{input_file_path}' "
                                f"does not exist.")

    if input_file_path.lower().endswith('.json'):
        with open(input_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            folder_paths = data.get('folder_paths', [])
    else:
        raise ValueError("Input file must be in .json format.")

    if not folder_paths:
        raise ValueError("Input file does not contain folder paths.")

    # Check if folders exist and print information about files
    valid_folder_paths = []
    for folder_path in folder_paths:
        if os.path.isdir(folder_path):
            files = os.listdir(folder_path)
            num_files = len(files)
            file_types = set([os.path.splitext(f)[1].lower() for f in files])
            print(f"\nFolder: {folder_path}")
            print(f"Number of files: {num_files}")
            print(f"File types: {', '.join(file_types)}")
            valid_folder_paths.append(folder_path)
        else:
            logging.warning(f"Folder '{folder_path}' does not exist.")

    if not valid_folder_paths:
        raise ValueError("No available folders for processing.")

    print(f"\nFound {len(valid_folder_paths)} available "
                 f"folders for processing.")
    return valid_folder_paths


def create_results_folders(folder_path, angle_value_str, timestamp):
    """
    Creates necessary folders to save results in the specified folder.

    Args:
        folder_path (str): Path to the input folder.
        angle_value_str (str): String representation
        of the angle (e.g., '15' or '15_5').
        timestamp (str): Current timestamp for naming folders.

    Returns:
        Tuple[str, str, str]: Paths to results, tables, and images folders.
    """
    results_folder_name = (f"Alignment_assay_results_angle"
                           f"_{angle_value_str}_{timestamp}")
    results_folder = os.path.join(folder_path, results_folder_name)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    logging.info(f"Results will be saved in: {results_folder}")

    # Add logging
    log_file = os.path.join(results_folder, 'log.log')
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logging.getLogger('').addHandler(file_handler)

    table_folder = os.path.join(results_folder, "Tables")
    images_folder = os.path.join(results_folder, "Images")
    Path(table_folder).mkdir(parents=True, exist_ok=True)
    Path(images_folder).mkdir(parents=True, exist_ok=True)
    print(f"Tables and Images folders created in {results_folder}")

    return results_folder, table_folder, images_folder


def process_part1(folder_path,
                  results_folder,
                  fibronectin_channel_index,
                  desired_width,
                  desired_height,
                  ij):
    """
    Part 1: Create 2D projections for images in the folder.

    Args:
        folder_path (str): Path to the input folder.
        results_folder (str): Path to the result folder.
        fibronectin_channel_index (int): Index of the
        fibronectin channel (starting from 1).
        desired_width (int): Desired width of output images.
        desired_height (int): Desired height of output images.
        ij: ImageJ context.

    Returns:
        Dict[str, Dict]: Information about Z-stacks processed for each folder.
    """
    print("\nPart 1: Processing images...")
    z_stacks_info_folder = {}

    # Import IJ and ZProjector
    IJ = sj.jimport('ij.IJ')
    ZProjector = sj.jimport('ij.plugin.ZProjector')

    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(('.tif', '.tiff', '.nd2')):
            logging.warning(f"Skipping file '{filename}' as it is not .nd2 or .tif.")
            continue

        file_path = os.path.join(folder_path, filename)
        print(f"\nProcessing file: {file_path}")

        # Close all windows before starting processing
        IJ.run("Close All")

        # Open the image using Bio-Formats
        imp = IJ.openImage(file_path)
        if imp is None:
            logging.warning(f"Failed to open image '{file_path}'. Check Bio-Formats plugin.")
            continue

        # Get image dimensions
        width, height, channels, slices, frames = imp.getDimensions()
        print(f"Image dimensions for '{filename}': "
                     f"width={width}, "
                     f"height={height}, "
                     f"channels={channels}, "
                     f"slices={slices}, "
                     f"frames={frames}")

        # Check if the specified channel is available
        if fibronectin_channel_index > channels:
            print(f"Specified channel exceeds "
                         f"available channels in '{filename}'. "
                         f"Skipping file.")
            imp.close()
            continue

        # Process the fibronectin channel
        print(f"Processing fibronectin channel "
                     f"({fibronectin_channel_index}) in '{filename}'.")
        imp.setC(fibronectin_channel_index)
        IJ.run(imp, "Duplicate...",
               f"title=imp_fibro duplicate "
               f"channels={fibronectin_channel_index}")
        imp_fibro = IJ.getImage()

        # Perform maximum intensity projection along Z
        zp_fibro = ZProjector(imp_fibro)
        zp_fibro.setMethod(ZProjector.MAX_METHOD)
        zp_fibro.doProjection()
        fibro_proj = zp_fibro.getProjection()
        fibro_proj = fibro_proj.resize(desired_width,
                                       desired_height,
                                       "bilinear")
        IJ.run(fibro_proj, "8-bit", "")  # Convert to grayscale

        output_filename = os.path.splitext(filename)[0] + '_processed.tif'
        output_path = os.path.join(results_folder, output_filename)
        IJ.saveAs(fibro_proj, "Tiff", output_path)
        logging.info(f"Processed image saved at '{output_path}'.")
        fibro_proj.close()
        imp_fibro.close()

        # Close original image
        imp.close()

        # Close all windows
        IJ.run("Close All")

        # Save Z-stack information
        processed_base_name = os.path.splitext(output_filename)[0]
        z_stacks_info_folder[processed_base_name] = {
            'original_filename': filename,
            'number_of_z_stacks': slices,
            'z_stack_type': 'slices'
        }

    print(f"\nPart 1 successfully completed for folder {folder_path}.")
    return z_stacks_info_folder


def process_part2_orientationpy(results_folder, images_folder):
    """
    Part 2: Applying orientationpy to 2D
    projections of the fibronectin channel.
    """
    print("\nPart 2: Applying orientationpy to "
                 "2D projections of the fibronectin channel....")
    processed_files = [f for f in os.listdir(results_folder)
                       if f.lower().endswith('_processed.tif')]

    if not processed_files:
        logging.warning("Processed images not found. "
                        "Make sure Part 1 was successfully completed.")
        return

    for filename in processed_files:
        file_path = os.path.join(results_folder, filename)
        print(f"\nProcessing file: {file_path}")

        # Read the image into a NumPy array and convert to float
        image_gray = io.imread(file_path).astype(float)
        print(f"Image '{filename}' successfully "
                     f"read with dimensions "
                     f"{image_gray.shape} in grayscale, "
                     f"with max value: {image_gray.max()}.")
    
        # Set correct anisotropy values
        anisotropy = np.array([1., 1., 1.])  # Relative pixel size for all axes
        gradient_mode = 'splines'
        gradients = computeGradient(image_gray,
                                    mode=gradient_mode,
                                    anisotropy=anisotropy)
        # Gradients computed
        print(f"Gradients for '{filename}' "
                     f"computed using mode {gradient_mode} "
                     f"and anisotropy {anisotropy}.")

        # # Visualisation of Gx и Gy gradients optional
        # plt.figure(figsize=(10, 4))
        # Gy, Gx = gradients
        # plt.subplot(1, 2, 1)
        # plt.title("Gy - Gradient Y")
        # plt.imshow(Gy, cmap="coolwarm", vmin=-64, vmax=64)
        # plt.colorbar(shrink=0.7)

        # plt.subplot(1, 2, 2)
        # plt.title("Gx - Gradient X")
        # plt.imshow(Gx, cmap="coolwarm", vmin=-64, vmax=64)
        # plt.colorbar(shrink=0.7)

        # plt.tight_layout()
        # gradient_path = os.path.join(images_folder, f"{os.path.splitext(filename)[0]}_gradients.png")
        # plt.savefig(gradient_path)
        # plt.close()
        # print(f"Gradients saved in '{gradient_path}'.")

        # Compute the structure tensor with appropriate sigma
        sigma = 1  # Reduce sigma if needed
        structure_tensor = computeStructureTensor(gradients, sigma=sigma)
        directionality = computeStructureDirectionality(structure_tensor)
        orientations = computeOrientation(structure_tensor)

        print(f"Structure Tensor, "
                     f"intensity, "
                     f"directionality, "
                     f"and orientation computed "
                     f"for '{filename}' with "
                     f"sigma parameter {sigma}.")

        # # Check the value of direction
        # print(f"Values of direction: min {directionality.min()},bmax {directionality.max()}, mean {directionality.mean()}, st.dev {directionality.std()}")

        # Normalize directionality for visualization
        vmin, vmax = 10, 1e8
        normalized_directionality = np.clip(directionality, vmin, vmax)
        normalized_directionality = np.log(normalized_directionality)
        normalized_directionality -= normalized_directionality.min()
        normalized_directionality /= normalized_directionality.max()
        normalized_directionality[image_gray == 0] = 0

        # # Visualization of intensity and directionality (optional)
        # plt.figure(figsize=(10, 4))

        # plt.subplot(1, 2, 1)
        # plt.imshow(intensity / intensity.max(), vmin=0, vmax=1)
        # plt.colorbar(shrink=0.7)
        # plt.title("Intensity Normalized")

        # plt.subplot(1, 2, 2)
        # plt.imshow(directionality, norm=matplotlib.colors.LogNorm(vmin=10, vmax=1e8))
        # plt.colorbar(shrink=0.7)
        # plt.title("Directionality Normalized")

        # plt.tight_layout()
        # intensity_directionality_path = os.path.join(images_folder, f"{os.path.splitext(filename)[0]}_intensity_directionality.png")
        # plt.savefig(intensity_directionality_path)
        # plt.close()
        # print(f"Intensity and directionality images saved at '{intensity_directionality_path}'.")

        # # Visualization of overlaying orientations on the image
        # plt.figure(figsize=(6, 6))
        # plt.imshow(image_gray, cmap="Greys_r", vmin=0)

        # plt.imshow(orientations["theta"], cmap="hsv", alpha=0.5, vmin=-90, vmax=90)

        # plt.colorbar(label='Angle (degrees)')
        # plt.title(f"Overlay with Orientation for {filename}")
        # overlay_path = os.path.join(images_folder, f"{os.path.splitext(filename)[0]}_orientation_overlay.png")
        # plt.savefig(overlay_path)
        # plt.close()
        # print(f"Orientation overlay image saved at '{overlay_path}'.")

        # Optional analysis using box pixels. If used, replace in orientation_flat with orientations to orientationsBoxes
        # boxSizePixels = 3
        # structureTensorBoxes = computeStructureTensorBoxes(
        #     [Gy, Gx],
        #     [boxSizePixels, boxSizePixels],
        # )

        # orientationsBoxes = computeOrientation(
        #     structureTensorBoxes,
        #     mode="fiber",
        # )

        # Create histogram of orientation distribution
        orientation_flat = orientations["theta"].flatten()
        orientation_flat = orientation_flat[~np.isnan(orientation_flat)]

        print(f"Calculating orientation histogram for '{filename}'.")
        num_bins = 180  # For 1-degree bins
        hist, bin_edges = np.histogram(orientation_flat, bins=num_bins, range=(-90, 90))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Prepare DataFrame for saving
        df = pd.DataFrame({
            'orientation_angle': bin_centers,
            'occurrence_value': hist
        })

        # Save orientation distribution to CSV
        table_folder = os.path.join(results_folder, "Tables")
        os.makedirs(table_folder, exist_ok=True)
        excel_filename = os.path.splitext(filename)[0] + '_orientation_distribution.csv'
        excel_path = os.path.join(table_folder, excel_filename)
        df.to_csv(excel_path, index=False)
        logging.info(f"Orientation distribution data saved at '{excel_path}'.")

        # Generate orientation composition image (HSV composition)
        print(f"Generating orientation composition image for '{filename}'.")

        imDisplayHSV = np.zeros((image_gray.shape[0], image_gray.shape[1], 3), dtype="f4")
        # Hue = (angle + 90)/180 -> normalized to [0, 1]
        imDisplayHSV[:, :, 0] = (orientations["theta"] + 90) / 180.0  
        # Saturation = normalized directionality
        imDisplayHSV[:, :, 1] = normalized_directionality
        # Value = original image normalized
        imDisplayHSV[:, :, 2] = image_gray / image_gray.max()

        # Create a figure and axes for the image
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(matplotlib.colors.hsv_to_rgb(imDisplayHSV))
        ax.set_title(f"Image-Orientation Composition for\n{filename}")

        # Create a ScalarMappable for the colorbar to represent angles
        sm = (matplotlib.cm
              .ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-90,
                                                               vmax=90),
                              cmap='hsv'))
        sm.set_array([])  # Required for proper colorbar functionality

        # Attach the colorbar to the figure and axes
        fig.colorbar(sm,
                     ax=ax,
                     orientation="vertical",
                     label="Degrees from Horizontal",
                     shrink=0.7)

        composition_path = os.path.join(images_folder,
                                        f"{os.path.splitext(filename)[0]}"
                                        f"_orientation_composition.png")
        fig.savefig(composition_path)
        plt.close(fig)
        logging.info(f"Orientation composition image saved at "
                     f"'{composition_path}'.")

        # Saving directionality values for debugging (optional)
        # directionality_path = os.path.join(images_folder, f"{os.path.splitext(filename)[0]}_directionality.npy")
        # np.save(directionality_path, normalized_directionality)
        # print(f"Normalized directionality values saved for debugging at '{directionality_path}'.")


def process_part3(results_folder,
                  analysis_folder,
                  angle_value,
                  z_stacks_info):
    """
    Part 3: Processing CSV files and summarizing results.

    Arguments:
        results_folder (str): Path to the result folder.
        analysis_folder (str): Path to the analysis folder.
        angle_value (float): Angle for analysis.
        z_stacks_info (Dict[str, Dict]):
        Information about Z-stacks processed for each folder.
    """
    print("\nPart 3: Processing CSV files "
                 "and generating summary of results...")

    table_folder = os.path.join(results_folder, "Tables")
    if not os.path.exists(table_folder):
        logging.warning(f"Folder '{table_folder}' does not exist. "
                        f"Make sure Part 2 was completed successfully.")
        return

    # Get list of CSV files in the Tables folder
    file_list = [f for f in os.listdir(table_folder) if f.endswith('.csv')]

    if not file_list:
        logging.warning(f"No CSV files found in '{table_folder}'. "
                        f"Make sure Part 2 was completed successfully.")
        return

    # Initialize list to store summary data
    summary_data = []

    # Process each CSV file
    for file_name in file_list:
        file_path = os.path.join(table_folder, file_name)
        print(f"\nProcessing CSV file: {file_name}")

        # Read CSV file into DataFrame
        read_file = pd.read_csv(file_path)

        # Rename columns
        read_file.rename(columns={read_file.columns[0]: "orientation_angle",
                                  read_file.columns[1]: "occurrence_value"},
                         inplace=True)

        # Find the corresponding angle
        angle_of_max_occurrence_value = (read_file['orientation_angle'][read_file['occurrence_value'].idxmax()])

        # Normalize angles relative to the angle of maximum value
        read_file['angles_normalized_to_angle_of_MOV'] = (read_file['orientation_angle']
                                                          - angle_of_max_occurrence_value)

        read_file['corrected_angles'] = (read_file['angles_normalized_to_angle_of_MOV']
                                         .apply(correct_angle))

        # Rank corrected angles
        read_file['rank_of_angle_occurrence_value'] = read_file['corrected_angles'].rank(method='min')

        # Sum of "occurrence_value" values
        sum_of_occurrence_values = read_file['occurrence_value'].sum()

        # Calculate percentage
        read_file['percent_occurrence_value_to_sum_of_occurrence_value'] = (read_file['occurrence_value']
                                                                            / sum_of_occurrence_values) * 100

        # Filter rows by angle range
        filtered_data = read_file[(read_file['corrected_angles'] >= -angle_value) & (read_file['corrected_angles'] <= angle_value)]

        # Sum percentage values for filtered rows
        percentage_of_fibers_aligned_within_angle = filtered_data['percent_occurrence_value_to_sum_of_occurrence_value'].sum()

        # Determine orientation mode
        orientation_mode = 'disorganized' if percentage_of_fibers_aligned_within_angle < 55 else 'aligned'

        # Sort DataFrame by rank
        read_file_sorted = read_file.sort_values(by='rank_of_angle_occurrence_value')

        # Generate output file name and save sorted DataFrame
        output_file_name = f'{os.path.splitext(file_name)[0]}_processed.csv'
        output_file_path = os.path.join(analysis_folder, output_file_name)
        read_file_sorted.to_csv(output_file_path, index=False)
        logging.info(f"Processed data saved at: {output_file_path}")

        # Get the base name of the processed file
        processed_base_name = os.path.splitext(file_name.replace('_orientation_distribution.csv', ''))[0]

        # Get Z-stack information
        z_stacks_info_entry = z_stacks_info.get(processed_base_name, None)
        if z_stacks_info_entry is not None:
            number_of_z_stacks = z_stacks_info_entry['number_of_z_stacks']
            z_stack_type = z_stacks_info_entry['z_stack_type']
        else:
            number_of_z_stacks = 'N/A'
            z_stack_type = 'N/A'

        # Add data to summary list
        summary_data.append({
            'File_Name': file_name,
            'Number_of_Z_Stacks': number_of_z_stacks,
            'Z_Stack_Type': z_stack_type,
            f'Percentage_Fibers_Aligned_Within_'
            f'{angle_value}_Degree': percentage_of_fibers_aligned_within_angle,
            'Orientation_Mode': orientation_mode
        })

    # Save summary data to CSV
    summary_df = pd.DataFrame(summary_data)
    summary_file_path = os.path.join(analysis_folder,
                                     f'Alignment_Summary.csv')
    summary_df.to_csv(summary_file_path, index=False)
    logging.info(f"\nSummary data saved at: {summary_file_path}")

    logging.info(f"\nProcessing completed for folder {results_folder}. "
                 f"All results saved in folder: {results_folder}")


def process_folder(folder_path,
                   fibronectin_channel_index,
                   angle_value,
                   desired_width,
                   desired_height,
                   ij):
    """
    Processes one folder: Part 1, Part 2, and Part 3.

    Arguments:
        folder_path (str): Path to the folder for processing.
        fibronectin_channel_index (int): Index of the
        fibronectin channel for this folder.
        angle_value (float): Angle for analysis.
        desired_width (int): Desired width of output images.
        desired_height (int): Desired height of output images.
        ij: ImageJ context.
    """
    if not os.path.exists(folder_path):
        logging.warning(f"Folder '{folder_path}' does not exist. Skipping this folder.")
        return

    # Convert angle value to string for folder naming
    angle_str = f"{angle_value}".replace('.', '_')

    # Get current date and time for folder naming
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")

    # Create folders for results
    results_folder, table_folder, images_folder = create_results_folders(folder_path,
                                                                         angle_str,
                                                                         timestamp)

    # Part 1: Create 2D projections
    z_stacks_info_folder = process_part1(folder_path,
                                         results_folder,
                                         fibronectin_channel_index,
                                         desired_width,
                                         desired_height,
                                         ij)

    # Part 2: Orientation analysis using orientationpy
    process_part2_orientationpy(results_folder,
                                images_folder)

    # Part 3: Process CSV files and generate summary of results
    analysis_folder = os.path.join(results_folder, 'Analysis')
    if not os.path.exists(analysis_folder):
        os.makedirs(analysis_folder)
    process_part3(results_folder,
                  analysis_folder,
                  angle_value,
                  z_stacks_info_folder)


def main_fibronectin_processing(input_file_path,
                                angle_value=15,
                                desired_width=500,
                                desired_height=500):
    """
    Main function: Parse command line arguments and start processing.
    """
    # Setting up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Initialize ImageJ
    print("Initializing ImageJ in headless mode...")
    try:
        # Use mode='headless' instead of headless=True
        ij = imagej.init('sc.fiji:fiji', mode='headless')
    except Exception as e:
        raise ImageJInitializationError(f"Failed to initialize ImageJ: {e}")
    print("ImageJ successfully initialized.")

    # Get folder paths
    folder_paths = get_folder_paths(input_file_path)

    # Get fibronectin channel index
    fibronectin_channel_index = int(input("Enter fibronectin "
                                          "channel index (starting from 1): ").strip())

    # Prompt user to start processing
    start_processing = (input("\nDo you want to start processing the files? (yes/no): ")
                        .strip()
                        .lower())
    if start_processing not in ('yes', 'y'):
        raise ValueError("File processing canceled by user.")

    # Process each folder
    for folder_path in folder_paths:
        process_folder(folder_path,
                       fibronectin_channel_index,
                       angle_value,
                       desired_width,
                       desired_height,
                       ij)

    print("\nAll folders have been processed.")

    # Terminate ImageJ
    print("Terminating ImageJ...")
    ij.context().dispose()
    print("Script execution completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     'Orientation analysis script')
    parser.add_argument('-i',
                        '--input',
                        required=True,
                        help='Path to the input_paths.json '
                             'file containing paths to folders for processing')
    parser.add_argument('-a',
                        '--angle_value',
                        default=15,
                        help='Angle value')

    args = parser.parse_args()

    main_fibronectin_processing(args.input, args.angle_value)
