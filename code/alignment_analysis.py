#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import json
import logging
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
from orientationpy import (computeGradient, computeOrientation,
                           computeStructureDirectionality,
                           computeStructureTensor)
from skimage import io


class ImageJInitializationError(Exception):
    """
    Exception raised for unsuccessful initialization of ImageJ.
    """
    pass


def initialize_imagej():
    """
    Initialize ImageJ in headless mode.

    Returns:
        ij (imagej.ImageJ): The initialized ImageJ instance.
    """
    # Attempt to initialize ImageJ headless mode
    print("Initializing ImageJ...")
    try:
        ij = imagej.init('sc.fiji:fiji', mode='headless')
    except Exception as e:
        raise ImageJInitializationError(
            f"Failed to initialize ImageJ: {e}")
    print("ImageJ initialization completed.")
    return ij


def correct_angle(x):
    """
    Correct angles to the range of -90 to 90 degrees.

    Args:
        x (float): The angle value to correct.

    Returns:
        float: The corrected angle.
    """
    if x < -90:
        return x + 180
    elif x > 90:
        return x - 180
    else:
        return x


def get_folder_paths(input_file_path):
    """
    Read an input JSON file containing folder paths.

    Args:
        input_file_path (str): Path to the input JSON file.

    Returns:
        List[str]: List of valid folder paths.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file does not contain folder paths or
                    no valid folders exist.
    """
    if not os.path.isfile(input_file_path):
        raise FileNotFoundError(
            f"File '{input_file_path}' does not exist."
        )

    if input_file_path.lower().endswith('.json'):
        with open(input_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
            folder_paths = data.get('folder_paths', [])
    else:
        raise ValueError("Input file must be in .json format.")

    if not folder_paths:
        raise ValueError("Input file does not contain folder paths.")

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
            raise ValueError(f"Folder '{folder_path}' does not exist.")

    if not valid_folder_paths:
        raise ValueError("No available folders for processing.")

    print(f"\nFound {len(valid_folder_paths)} "
          f"available folders for processing.")
    return valid_folder_paths


def create_results_folders(folder_path, angle_value_str, timestamp):
    """
    Create folders to save results in the specified folder.

    Args:
        folder_path (str): Path to the input folder.
        angle_value_str (str): String representation of the angle.
        timestamp (str): Current timestamp for naming folders.

    Returns:
        Tuple[str, str, str]: Paths to results, tables, and images folders.
    """
    results_folder_name = (
        f"Alignment_assay_results_angle_{angle_value_str}_{timestamp}"
    )
    results_folder = os.path.join(folder_path, results_folder_name)
    Path(results_folder).mkdir(parents=True, exist_ok=True)

    # Add logging
    file_handler = logging.FileHandler(os.path.join(results_folder,
                                                    'log.log'),
                                       mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(
        logging.Formatter('%(asctime)s '
                          '- %(levelname)s '
                          '- %(message)s'))
    logging.getLogger('').addHandler(file_handler)

    logging.info(f"Results will be saved in: {results_folder}")

    table_folder = os.path.join(results_folder, "Tables")
    images_folder = os.path.join(results_folder, "Images")
    Path(table_folder).mkdir(parents=True, exist_ok=True)
    Path(images_folder).mkdir(parents=True, exist_ok=True)
    logging.info(f"Tables and Images folders "
                 f"created in {results_folder}")

    return results_folder, table_folder, images_folder


def process_part1(
    folder_path,
    results_folder,
    fibronectin_channel_index,
    desired_width,
    desired_height,
    ij
):
    """
    Part 1: Create 2D projections for images in the folder.

    Args:
        folder_path (str): Path to the input folder.
        results_folder (str): Path to the result folder.
        fibronectin_channel_index (int): Index of the fibronectin channel.
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
        if (
            filename.startswith("._") or
            filename.startswith(".") or
            not filename.lower().endswith(('.tif', '.tiff', '.nd2'))
        ):
            logging.warning(f"Skipping hidden or invalid file: '{filename}'")
            continue

        file_path = os.path.join(folder_path, filename)
        logging.info(f"\nProcessing file: {file_path}")

        # Close all windows before starting processing
        IJ.run("Close All")

        # Open the image using Bio-Formats
        imp = IJ.openImage(file_path)
        if imp is None:
            logging.warning(
                f"Failed to open image '{file_path}'. "
                f"Check Bio-Formats plugin."
            )
            continue

        # Get image dimensions
        width, height, channels, slices, frames = imp.getDimensions()
        logging.info(
            f"Image dimensions for '{filename}': width={width}, "
            f"height={height}, channels={channels}, slices={slices}, "
            f"frames={frames}"
        )

        # Check if the specified channel is available
        if fibronectin_channel_index > channels:
            logging.warning(
                f"Specified channel exceeds "
                f"available channels in '{filename}'. "
                f"Skipping file."
            )
            imp.close()
            continue

        # Process the fibronectin channel
        logging.info(
            f"Processing fibronectin "
            f"channel ({fibronectin_channel_index}) "
            f"in '{filename}'."
        )
        imp.setC(fibronectin_channel_index)
        IJ.run(
            imp,
            "Duplicate...",
            f"title=imp_fibro duplicate channels={fibronectin_channel_index}"
        )
        imp_fibro = IJ.getImage()

        # Perform maximum intensity projection along Z
        zp_fibro = ZProjector(imp_fibro)
        zp_fibro.setMethod(ZProjector.MAX_METHOD)
        zp_fibro.doProjection()
        fibro_proj = zp_fibro.getProjection()
        fibro_proj = fibro_proj.resize(desired_width, desired_height,
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


def process_part2_orientationpy(results_folder,
                                images_folder):
    """
    Part 2: Apply orientationpy to 2D projections of the fibronectin channel.

    Args:
        results_folder (str): Path to the folder with processed images.
        images_folder (str): Path to the folder where images will be saved.
    """
    print(
        "\nPart 2: Applying orientationpy to 2D projections "
        "of the fibronectin channel..."
    )

    # Create subfolder for normalized images
    normalized_images_folder = os.path.join(images_folder,
                                            "normalized_images")
    Path(normalized_images_folder).mkdir(parents=True,
                                         exist_ok=True)
    processed_files = [
        f for f in os.listdir(results_folder)
        if (f.lower().endswith('_processed.tif')
            and not f.startswith("._")
            and not f.startswith("."))
    ]

    if len(processed_files) == 0:
        logging.warning(
            "Processed images not found. "
            "Make sure Part 1 was completed."
        )
        return

    for filename in processed_files:
        file_path = os.path.join(results_folder, filename)
        logging.info(f"\nProcessing file: {file_path}")

        # Read the image into a NumPy array and convert to float
        image_gray = io.imread(file_path).astype(float)
        logging.info(
            f"Image '{filename}' successfully read with dimensions "
            f"{image_gray.shape}, max value: {image_gray.max()}."
        )

        # Compute gradients
        anisotropy = np.array([1., 1., 1.])  # Relative pixel size
        gradient_mode = 'splines'
        gradients = computeGradient(image_gray,
                                    mode=gradient_mode,
                                    anisotropy=anisotropy)
        logging.info(
            f"Gradients for '{filename}' computed "
            f"using mode {gradient_mode}, anisotropy {anisotropy}."
        )

        # Compute the structure tensor
        sigma = 3  # Standard deviation for Gaussian
        structure_tensor = computeStructureTensor(gradients,
                                                  sigma=sigma)
        directionality = computeStructureDirectionality(structure_tensor)
        orientations = computeOrientation(structure_tensor)

        logging.info(
            f"Structure Tensor, intensity, directionality, and orientation "
            f"computed for '{filename}' with sigma={sigma}."
        )

        # Normalize directionality
        vmin, vmax = 10, 1e8
        normalized_directionality = np.clip(directionality,
                                            vmin,
                                            vmax)
        normalized_directionality = np.log(normalized_directionality)
        normalized_directionality -= normalized_directionality.min()
        normalized_directionality /= normalized_directionality.max()
        normalized_directionality[image_gray == 0] = 0

        # Create histogram of orientation distribution
        # orientation_flat = orientations["theta"].flatten()
        # orientation_flat = orientation_flat[~np.isnan(orientation_flat)]

        # print(f"Calculating orientation histogram for '{filename}'.")
        # num_bins = 180  # 1-degree bins
        # hist, bin_edges = np.histogram(orientation_flat,
        # bins=num_bins, range=(-90, 90))
        # bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Calculate the histogram in OrientationJ-compatible way
        theta = orientations["theta"]  # grad
        gx, gy = np.gradient(image_gray)  # gradients
        energy = gx**2 + gy**2  # Weight = Energy

        mask = energy > 0  # as OrientationJ
        theta_flat = theta[mask].ravel()
        energy_flat = energy[mask].ravel()

        num_bins = 180  # 1°-bin
        hist, bin_edges = np.histogram(
            theta_flat,
            bins=num_bins,
            range=(-90, 90),
            weights=energy_flat
        )

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        df = pd.DataFrame({
            "ori_angle": bin_centers,
            "occ_value": hist
        })

        # Save orientation distribution to CSV
        table_folder = os.path.join(results_folder, "Tables")
        os.makedirs(table_folder, exist_ok=True)
        excel_filename = (f"{os.path.splitext(filename)[0]}"
                          f"_orientation_distribution.csv")
        excel_path = os.path.join(table_folder, excel_filename)
        df.to_csv(excel_path, index=False)
        logging.info(f"Orientation distribution data saved at '{excel_path}'.")

        # Generate orientation composition image (HSV)
        logging.info(f"Generating orientation "
                     f"composition image for '{filename}'.")
        im_display_hsv = np.zeros((image_gray.shape[0],
                                   image_gray.shape[1], 3),
                                  dtype="f4")

        # Hue: (angle + 90)/180 -> normalized to [0, 1]
        im_display_hsv[:, :, 0] = (orientations["theta"] + 90) / 180.0
        # Saturation: normalized directionality
        im_display_hsv[:, :, 1] = normalized_directionality
        # Value: original image normalized
        im_display_hsv[:, :, 2] = image_gray / image_gray.max()

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(matplotlib.colors.hsv_to_rgb(im_display_hsv))
        ax.set_title(f"Image-Orientation Composition for\n{filename}")

        # Create ScalarMappable for the colorbar to represent angles
        sm = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=-90, vmax=90),
            cmap='hsv'
        )
        sm.set_array([])

        fig.colorbar(
            sm,
            ax=ax,
            orientation="vertical",
            label="Degrees from Horizontal",
            shrink=0.7
        )

        composition_path = os.path.join(
            images_folder,
            f"{os.path.splitext(filename)[0]}_orientation_composition.png"
        )
        fig.savefig(composition_path)
        plt.close(fig)
        logging.info(f"Orientation composition image "
                     f"saved at '{composition_path}'.")

        # Calculate modal angle (dominant orientation)
        modal_angle = bin_centers[np.argmax(hist)]
        logging.info(f"Modal angle for {filename}: {modal_angle:.2f}°")

        # Create HSV representation
        im_display_hsv = np.zeros((image_gray.shape[0],
                                   image_gray.shape[1], 3),
                                  dtype="f4")
        im_display_hsv[:, :, 0] = (orientations["theta"] + 90) / 180.0
        im_display_hsv[:, :, 1] = normalized_directionality
        im_display_hsv[:, :, 2] = image_gray / image_gray.max()

        # Normalize Hue: shift to display modal angle as cyan (180°)
        hue_shift = -2 * modal_angle  # Required shift in degrees
        normalized_hue = (im_display_hsv[:, :, 0] + hue_shift/360) % 1.0
        im_display_hsv_normalized = im_display_hsv.copy()
        im_display_hsv_normalized[:, :, 0] = normalized_hue

        # Save original orientation composition
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(matplotlib.colors.hsv_to_rgb(im_display_hsv))
        ax.set_title(f"Image-Orientation Composition for\n{filename}")
        sm = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=-90, vmax=90),
            cmap='hsv'
        )
        sm.set_array([])
        fig.colorbar(
            sm,
            ax=ax,
            orientation="vertical",
            label="Degrees from Horizontal",
            shrink=0.7
        )
        composition_path = os.path.join(
            images_folder,
            f"{os.path.splitext(filename)[0]}_orientation_composition.png"
        )
        fig.savefig(composition_path)
        plt.close(fig)

        # Save normalized orientation composition
        fig_norm, ax_norm = plt.subplots(figsize=(6, 6))
        ax_norm.imshow(matplotlib.colors.hsv_to_rgb(im_display_hsv_normalized))
        ax_norm.set_title(f"Normalized Orientation for\n{filename}")
        sm_norm = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(vmin=-90, vmax=90),
            cmap='hsv'
        )
        sm_norm.set_array([])
        fig_norm.colorbar(
            sm_norm,
            ax=ax_norm,
            orientation="vertical",
            label="Deviation from Dominant Direction (°)",
            shrink=0.7
        )
        normalized_path = os.path.join(
            normalized_images_folder,
            f"{os.path.splitext(filename)[0]}_normalized_orientation.png"
        )
        fig_norm.savefig(normalized_path)
        plt.close(fig_norm)
        logging.info(f"Normalized orientation image "
                     f"saved at '{normalized_path}'.")


def process_part3(results_folder,
                  analysis_folder,
                  angle_value,
                  z_stacks_info):
    """
    Part 3: Process CSV files and summarize results.

    Args:
        results_folder (str): Path to the folder with results.
        analysis_folder (str): Path to the folder for analysis outputs.
        angle_value (float): Angle (in degrees) for analysis.
        z_stacks_info (Dict[str, Dict]): Z-stack information from Part 1.
    """
    print("\nPart 3: Processing CSV files and "
          "generating summary of results...")

    table_folder = os.path.join(results_folder, "Tables")
    if not os.path.exists(table_folder):
        logging.warning(
            f"Folder '{table_folder}' does not exist. "
            f"Make sure Part 2 was completed successfully."
        )
        return

    file_list = [
        f for f in os.listdir(table_folder)
        if (f.endswith('.csv')
            and not f.startswith("._")
            and not f.startswith("."))
    ]
    if not file_list:
        logging.warning(
            f"No CSV files found in '{table_folder}'. "
            f"Make sure Part 2 was completed successfully."
        )
        return

    summary_data = []

    for file_name in file_list:
        file_path = os.path.join(table_folder, file_name)
        logging.info(f"\nProcessing CSV file: {file_name}")

        # Read CSV file into DataFrame
        read_file = pd.read_csv(file_path)

        # Rename columns
        read_file.rename(
            columns={
                read_file.columns[0]: "ori_angle",
                read_file.columns[1]: "occ_value"
            },
            inplace=True
        )

        # Find the angle of maximum occupancy value
        angle_of_max_occ_value = read_file['ori_angle'][
            read_file['occ_value'].idxmax()
        ]

        # Normalize angles relative to the angle of maximum value
        read_file['angles_normalized_to_angle_of_MOV'] = (
            read_file['ori_angle'] - angle_of_max_occ_value
        )
        read_file['corrected_angles'] = read_file[
            'angles_normalized_to_angle_of_MOV'
        ].apply(correct_angle)

        # Rank corrected angles
        read_file['rank_of_angle_occ_value'] = read_file[
            'corrected_angles'
        ].rank(method='min')

        # Compute percentages
        sum_of_occ_values = read_file['occ_value'].sum()
        read_file['perc_occvalue2sum_of_occvalue'] = (
            read_file['occ_value'] / sum_of_occ_values
        ) * 100

        # Filter rows by angle range
        filtered_data = read_file[
            (read_file['corrected_angles'] >= -angle_value)
            & (read_file['corrected_angles'] <= angle_value)
        ]
        percentage_of_fibers_aligned_within_angle = (
            filtered_data['perc_occvalue2sum_of_occvalue'].sum()
        )

        # Determine orientation mode
        orientation_mode = 'disorganized'
        if percentage_of_fibers_aligned_within_angle >= 55:
            orientation_mode = 'aligned'

        # Sort DataFrame
        read_file_sorted = read_file.sort_values(by='rank_of_angle_occ_value')

        # Generate output file name
        output_file_name = f"{os.path.splitext(file_name)[0]}_processed.csv"
        output_file_path = os.path.join(analysis_folder, output_file_name)
        read_file_sorted.to_csv(output_file_path, index=False)
        logging.info(f"Processed data saved at: {output_file_path}")

        processed_base_name = os.path.splitext(
            file_name.replace('_orientation_distribution.csv', '')
        )[0]

        # Get Z-stack info
        z_stacks_info_entry = z_stacks_info.get(processed_base_name, None)
        if z_stacks_info_entry is not None:
            number_of_z_stacks = z_stacks_info_entry['number_of_z_stacks']
            z_stack_type = z_stacks_info_entry['z_stack_type']
        else:
            number_of_z_stacks = 'N/A'
            z_stack_type = 'N/A'

        summary_data.append({
            'File_Name': file_name,
            'Number_of_Z_Stacks': number_of_z_stacks,
            'Z_Stack_Type': z_stack_type,
            f'Percentage_Fibers_Aligned_Within_{angle_value}_Degree':
                percentage_of_fibers_aligned_within_angle,
            'Orientation_Mode': orientation_mode
        })

    # Save summary data
    summary_df = pd.DataFrame(summary_data)
    summary_file_path = os.path.join(analysis_folder,
                                     'Alignment_Summary.csv')
    summary_df.to_csv(summary_file_path, index=False)
    logging.info(f"\nSummary data saved at: {summary_file_path}")

    logging.info(
        f"\nProcessing completed for folder {results_folder}. "
        f"All results saved in folder: {results_folder}"
    )


def process_folder(
    folder_path,
    fibronectin_channel_index,
    angle_value,
    desired_width,
    desired_height,
    ij
):
    """
    Process a single folder (Part 1, Part 2, and Part 3).

    Args:
        folder_path (str): Path to the folder for processing.
        fibronectin_channel_index (int): Index of the fibronectin channel.
        angle_value (float): Angle (in degrees) for analysis.
        desired_width (int): Desired width of output images.
        desired_height (int): Desired height of output images.
        ij: ImageJ context.
    """
    if not os.path.exists(folder_path):
        print(
            f"Folder '{folder_path}' does not exist. Skipping this folder."
        )
        return

    # Convert angle value to string for folder naming
    angle_str = f"{angle_value}".replace('.', '_')

    # Create folders for results
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    results_folder, table_folder, images_folder = create_results_folders(
        folder_path,
        angle_str,
        timestamp
    )

    # Part 1: Create 2D projections
    z_stacks_info_folder = process_part1(
        folder_path,
        results_folder,
        fibronectin_channel_index,
        desired_width,
        desired_height,
        ij
    )

    # Part 2: Orientation analysis
    process_part2_orientationpy(results_folder, images_folder)

    # Part 3: Process CSV files and generate summary
    analysis_folder = os.path.join(results_folder, 'Analysis')
    if not os.path.exists(analysis_folder):
        os.makedirs(analysis_folder)

    process_part3(
        results_folder,
        analysis_folder,
        angle_value,
        z_stacks_info_folder
    )


def main_fibronectin_processing(
    input_file_path,
    angle_value=15,
    desired_width=500,
    desired_height=500
):
    """
    Main function to perform orientation-based
    analysis on fibronectin images.
    It:
    1) Reads folder paths from a JSON file.
    2) Creates 2D projections of fibronectin
    channels (Part 1).
    3) Applies orientation analysis using
    orientationpy (Part 2).
    4) Processes the resulting CSV files and
    generates a summary (Part 3).

    Args:
        input_file_path (str): Path to the JSON file with folder paths.
        angle_value (float): Angle for the alignment analysis.
        desired_width (int): Desired width of the 2D projected images.
        desired_height (int): Desired height of the 2D projected images.
    """
    logging.basicConfig(level=logging.INFO,
                        format='%(levelname)s: %(message)s')

    # Initialize ImageJ
    ij = initialize_imagej()

    # Get folder paths
    folder_paths = get_folder_paths(input_file_path)

    # Get fibronectin channel index
    fibr_chan_index = int(
        input("Enter fibronectin channel index (starting from 1): ").strip()
    )

    # Prompt user to start processing
    start_analysis = (input("\nDo you want to start processing? (y/n): ")
                      .strip().lower())
    if start_analysis in ('no', 'n'):
        ij.dispose()
        raise ValueError("Analysis canceled by user.")
    elif start_analysis not in ('yes', 'y', 'no', 'n'):
        raise ValueError("Incorrect input. Please enter y/n or yes/no")

    # Process each folder
    for folder_path in folder_paths:
        process_folder(
            folder_path,
            fibr_chan_index,
            angle_value,
            desired_width,
            desired_height,
            ij
        )

    print("\nAll folders have been processed.")
    print("Terminating ImageJ...")
    ij.context().dispose()
    print("Script execution completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Orientation analysis script')
    parser.add_argument(
        '-i',
        '--input',
        required=True,
        help='Path to the input_paths.json file containing '
             'paths to folders for processing'
    )
    parser.add_argument(
        '-a',
        '--angle_value',
        default=15,
        type=float,
        help='Angle value (in degrees)'
    )

    args = parser.parse_args()

    main_fibronectin_processing(args.input, args.angle_value)
