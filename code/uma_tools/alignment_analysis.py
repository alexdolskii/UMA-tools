#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import logging
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

from .config import read_config
from .files import assay_directory
from .image_run import ImageRun, active_run, image_attempt, image_names
from .imagej import (
    initialize_imagej,
)
from .plot_style import FONT_SIZES, PNG_DPI, rc_parameters, wrap_label
from .progress import (
    CANCELLATIONS,
    ask_channel,
    confirm_start,
    folder_error,
    folder_logged,
    folder_scope,
    outcome,
    phase,
    register_sources,
)
from .run import scoped_file_log, unique_output

# Select the backend before pyplot and image-processing imports.
matplotlib.use("Agg")
import matplotlib.colors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scyjava as sj  # noqa: E402
from orientationpy import (  # noqa: E402
    computeGradient,
    computeOrientation,
    computeStructureDirectionality,
    computeStructureTensor,
)
from skimage import io  # noqa: E402

PROGRESS_STAGES = ("Projection", "Orientation", "Summary")

_LOGGER = logging.getLogger(__name__)


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
    """Validate JSON; process missing folders as individual failures."""
    return [str(path) for path in read_config(Path(input_file_path))]


def create_results_folders(folder_path, angle_value_str, timestamp):
    """
    Create folders to save results in the specified folder.

    Args:
        folder_path (str): Path to the input folder.
        angle_value_str (str): String representation of the angle.
        timestamp (str): Current timestamp for naming folders.

    Returns:
        Tuple[str, str, str]: Paths to results, tables, and images
        folders.
    """
    _, output = unique_output(
        assay_directory(Path(folder_path), create=True),
        f"Alignment_assay_results_angle_{angle_value_str}_",
        timestamp=timestamp,
        include_pid=True,
    )
    results_folder = str(output)

    table_folder = os.path.join(results_folder, "Tables")
    images_folder = os.path.join(results_folder, "Images")
    Path(table_folder).mkdir(parents=True, exist_ok=True)
    Path(images_folder).mkdir(parents=True, exist_ok=True)
    _LOGGER.info(f"Tables and Images folders created in {results_folder}")

    return results_folder, table_folder, images_folder


def process_part1(
    folder_path,
    results_folder,
    fibronectin_channel_index,
    desired_width,
    desired_height,
    ij,
):
    """
    Part 1: Create 2D projections for images in the folder.

    Args:
        folder_path (str): Path to the input folder.
        results_folder (str): Path to the result folder.
        fibronectin_channel_index (int): Index of the fibronectin
        channel.
        desired_width (int): Desired width of output images.
        desired_height (int): Desired height of output images.
        ij: ImageJ context.

    Returns:
        Dict[str, Dict]: Information about Z-stacks processed for each
        folder.
    """
    phase(f"Stage 1/{len(PROGRESS_STAGES)}: {PROGRESS_STAGES[0]}")
    z_stacks_info_folder = {}

    # Import IJ and ZProjector
    IJ = sj.jimport("ij.IJ")
    ZProjector = sj.jimport("ij.plugin.ZProjector")
    Duplicator = sj.jimport("ij.plugin.Duplicator")

    for filename in image_names(folder_path, (".tif", ".tiff", ".nd2")):
        run = active_run()
        if run is not None and not run.eligible(filename):
            continue
        with image_attempt(filename, "Projection"):
            file_path = os.path.join(folder_path, filename)
            _LOGGER.info(f"\nProcessing file: {file_path}")

            imp = imp_fibro = fibro_proj = None
            try:
                # Close all windows before starting processing
                IJ.run("Close All")

                # Open the image using Bio-Formats
                imp = IJ.openImage(file_path)
                if imp is None:
                    raise ValueError(f"Could not open image: {file_path}")

                # Get image dimensions
                width, height, channels, slices, frames = imp.getDimensions()
                _LOGGER.info(
                    f"Image dimensions for '{filename}': width={width}, "
                    f"height={height}, channels={channels}, slices={slices}, "
                    f"frames={frames}"
                )

                # Check if the specified channel is available
                if not 1 <= fibronectin_channel_index <= channels:
                    imp.close()
                    raise ValueError(
                        f"Channel {fibronectin_channel_index} unavailable; "
                        f"image has {channels} channels"
                    )

                # Process the fibronectin channel
                _LOGGER.info(
                    f"Processing fibronectin "
                    f"channel ({fibronectin_channel_index}) "
                    f"in '{filename}'."
                )
                imp.setC(fibronectin_channel_index)
                imp_fibro = Duplicator().run(
                    imp,
                    fibronectin_channel_index,
                    fibronectin_channel_index,
                    1,
                    imp.getNSlices(),
                    1,
                    imp.getNFrames(),
                )
                imp_fibro.setTitle("imp_fibro")

                # Perform maximum intensity projection along Z
                zp_fibro = ZProjector(imp_fibro)
                zp_fibro.setMethod(ZProjector.MAX_METHOD)
                zp_fibro.doProjection()
                fibro_proj = zp_fibro.getProjection()
                fibro_proj = fibro_proj.resize(
                    desired_width, desired_height, "bilinear"
                )
                IJ.run(fibro_proj, "8-bit", "")  # Convert to grayscale

                output_filename = (
                    os.path.splitext(filename)[0] + "_processed.tif"
                )
                output_path = os.path.join(results_folder, output_filename)
                IJ.saveAs(fibro_proj, "Tiff", output_path)
                _LOGGER.info(f"Processed image saved at '{output_path}'.")
                fibro_proj.close()
                imp_fibro.close()

                # Close original image
                imp.close()

                # Close all windows
                IJ.run("Close All")

            finally:
                for image in (fibro_proj, imp_fibro, imp):
                    if image is not None:
                        image.close()
                IJ.run("Close All")

            # Save Z-stack information
            processed_base_name = os.path.splitext(output_filename)[0]
            z_stacks_info_folder[processed_base_name] = {
                "original_filename": filename,
                "number_of_z_stacks": slices,
                "z_stack_type": "slices",
            }

    phase(f"Stage 1/{len(PROGRESS_STAGES)}: Projection finished")
    return z_stacks_info_folder


@dataclass(frozen=True)
class OrientationResult:
    """
    Orientation arrays and the unchanged energy-weighted distribution.
    """

    orientations: dict
    normalized_directionality: np.ndarray
    bin_centers: np.ndarray
    histogram: np.ndarray
    distribution: pd.DataFrame


def calculate_orientation(
    image_gray: np.ndarray, filename: str
) -> OrientationResult:
    """
    Calculate the existing tensor, directionality, and 180-bin
    histogram.
    """
    # Compute gradients
    anisotropy = np.array([1.0, 1.0, 1.0])  # Relative pixel size
    gradient_mode = "splines"
    gradients = computeGradient(
        image_gray, mode=gradient_mode, anisotropy=anisotropy
    )
    _LOGGER.info(
        f"Gradients for '{filename}' computed "
        f"using mode {gradient_mode}, anisotropy {anisotropy}."
    )

    # Compute the structure tensor
    sigma = 3  # Standard deviation for Gaussian
    structure_tensor = computeStructureTensor(gradients, sigma=sigma)
    directionality = computeStructureDirectionality(structure_tensor)
    orientations = computeOrientation(structure_tensor)

    _LOGGER.info(
        f"Structure Tensor, intensity, directionality, and orientation "
        f"computed for '{filename}' with sigma={sigma}."
    )

    # Normalize directionality
    vmin, vmax = 10, 1e8
    normalized_directionality = np.clip(directionality, vmin, vmax)
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
        theta_flat, bins=num_bins, range=(-90, 90), weights=energy_flat
    )

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    df = pd.DataFrame({"ori_angle": bin_centers, "occ_value": hist})
    return OrientationResult(
        orientations, normalized_directionality, bin_centers, hist, df
    )


def save_orientation_result(
    result: OrientationResult,
    image_gray: np.ndarray,
    filename: str,
    results_folder: str,
    images_folder: str,
    normalized_images_folder: str,
    source_filename: str | None = None,
) -> None:
    """
    Export the original distribution and both orientation compositions.
    """
    orientations = result.orientations
    normalized_directionality = result.normalized_directionality
    bin_centers = result.bin_centers
    hist = result.histogram
    df = result.distribution
    # Save orientation distribution to CSV
    table_folder = os.path.join(results_folder, "Tables")
    os.makedirs(table_folder, exist_ok=True)
    excel_filename = (
        f"{os.path.splitext(filename)[0]}_orientation_distribution.csv"
    )
    excel_path = os.path.join(table_folder, excel_filename)
    df.to_csv(excel_path, index=False)
    _LOGGER.info(f"Orientation distribution data saved at '{excel_path}'.")

    # Generate orientation composition image (HSV)
    _LOGGER.info(f"Generating orientation composition image for '{filename}'.")
    im_display_hsv = np.zeros(
        (image_gray.shape[0], image_gray.shape[1], 3), dtype="f4"
    )

    # Hue: (angle + 90)/180 -> normalized to [0, 1]
    im_display_hsv[:, :, 0] = (orientations["theta"] + 90) / 180.0
    # Saturation: normalized directionality
    im_display_hsv[:, :, 1] = normalized_directionality
    # Value: original image normalized
    im_display_hsv[:, :, 2] = image_gray / image_gray.max()

    modal_angle = bin_centers[np.argmax(hist)]
    _LOGGER.info(f"Modal angle for {filename}: {modal_angle:.2f}°")
    # Preserve the original angular encoding and dominant-angle shift.
    hue_shift = -2 * modal_angle
    normalized_hsv = im_display_hsv.copy()
    normalized_hsv[:, :, 0] = (im_display_hsv[:, :, 0] + hue_shift / 360) % 1.0
    stem = os.path.splitext(filename)[0]
    for hsv, folder, suffix, title, label in (
        (
            im_display_hsv,
            images_folder,
            "orientation_composition",
            "Fiber orientation",
            "Orientation relative to horizontal (°)",
        ),
        (
            normalized_hsv,
            normalized_images_folder,
            "normalized_orientation",
            "Normalized fiber orientation",
            "Deviation from dominant orientation (°)",
        ),
    ):
        path = os.path.join(folder, f"{stem}_{suffix}.png")
        _save_orientation_figure(
            matplotlib.colors.hsv_to_rgb(hsv),
            source_filename or filename,
            path,
            title,
            label,
        )
        _LOGGER.info("Orientation composition saved at '%s'.", path)


def _save_orientation_figure(rgb, filename, path, title, colorbar_label):
    """Export the complete image with its title and literal filename."""
    width = 8.5
    footer = wrap_label(filename, (width - 0.6) * 72, FONT_SIZES["note"])
    footer_height = 0.24 * (footer.count("\n") + 1) + 0.12
    image_height = max(3.0, 6.6 * rgb.shape[0] / rgb.shape[1])
    with plt.rc_context(rc_parameters()):
        fig = plt.figure(
            figsize=(width, image_height + footer_height + 1),
            layout="constrained",
        )
        fig.get_layout_engine().set(w_pad=0.12, h_pad=0.12)
        try:
            grid = fig.add_gridspec(
                3, 1, height_ratios=[0.5, image_height, footer_height]
            )
            heading = fig.add_subplot(grid[0])
            heading.set_axis_off()
            heading.text(
                0.5,
                0.5,
                title,
                ha="center",
                va="center",
                fontsize=FONT_SIZES["title"],
                fontweight="bold",
            )
            body = grid[1].subgridspec(1, 2, width_ratios=[1, 0.045])
            axis = fig.add_subplot(body[0])
            axis.imshow(rgb, interpolation="nearest", aspect="equal")
            axis.set_axis_off()
            mapping = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(vmin=-90, vmax=90), cmap="hsv"
            )
            colorbar_axis = fig.add_subplot(body[1])
            colorbar = fig.colorbar(mapping, cax=colorbar_axis)
            colorbar.set_label(colorbar_label, fontsize=FONT_SIZES["axis"])
            colorbar.set_ticks([-90, -45, 0, 45, 90])
            colorbar.ax.tick_params(labelsize=FONT_SIZES["axis"])
            note = fig.add_subplot(grid[2])
            note.set_axis_off()
            note.text(
                0,
                0.5,
                footer,
                va="center",
                fontsize=FONT_SIZES["note"],
                transform=note.transAxes,
            )
            fig.savefig(path, dpi=PNG_DPI, facecolor="white")
        finally:
            plt.close(fig)


def process_part2_orientationpy(
    results_folder, images_folder, z_stacks_info=None
):
    """
    Part 2: Apply orientationpy to 2D projections of the fibronectin
    channel.

    Args:
        results_folder (str): Path to the folder with processed images.
        images_folder (str): Path to the folder where images will be
        saved.
    """
    phase(f"Stage 2/{len(PROGRESS_STAGES)}: {PROGRESS_STAGES[1]}")

    # Create subfolder for normalized images
    normalized_images_folder = os.path.join(images_folder, "normalized_images")
    Path(normalized_images_folder).mkdir(parents=True, exist_ok=True)
    processed_files = [
        f
        for f in os.listdir(results_folder)
        if (
            f.lower().endswith("_processed.tif")
            and not f.startswith("._")
            and not f.startswith(".")
            and os.path.isfile(os.path.join(results_folder, f))
        )
    ]

    if len(processed_files) == 0:
        _LOGGER.warning(
            "Processed images not found. Make sure Part 1 was completed."
        )
        return

    for filename in processed_files:
        name = (
            (z_stacks_info or {})
            .get(Path(filename).stem, {})
            .get("original_filename", filename)
        )
        run = active_run()
        if run is not None and not run.eligible(name):
            continue
        with image_attempt(name, "Orientation"):
            try:
                file_path = os.path.join(results_folder, filename)
                _LOGGER.info(f"\nProcessing file: {file_path}")

                # Read the image into a NumPy array and convert to float
                image_gray = io.imread(file_path).astype(float)
                _LOGGER.info(
                    f"Image '{filename}' successfully read with dimensions "
                    f"{image_gray.shape}, max value: {image_gray.max()}."
                )

                result = calculate_orientation(image_gray, filename)
                save_orientation_result(
                    result,
                    image_gray,
                    filename,
                    results_folder,
                    images_folder,
                    normalized_images_folder,
                    source_filename=name,
                )
            finally:
                plt.close("all")


def process_part3(results_folder, analysis_folder, angle_value, z_stacks_info):
    """
    Part 3: Process CSV files and summarize results.

    Args:
        results_folder (str): Path to the folder with results.
        analysis_folder (str): Path to the folder for analysis outputs.
        angle_value (float): Angle (in degrees) for analysis.
        z_stacks_info (Dict[str, Dict]): Z-stack information from Part
        1.
    """
    phase(f"Stage 3/{len(PROGRESS_STAGES)}: {PROGRESS_STAGES[2]}")

    table_folder = os.path.join(results_folder, "Tables")
    if not os.path.exists(table_folder):
        _LOGGER.warning(
            f"Folder '{table_folder}' does not exist. "
            f"Make sure Part 2 was completed successfully."
        )
        return

    file_list = [
        f
        for f in os.listdir(table_folder)
        if (
            f.endswith(".csv")
            and not f.startswith("._")
            and not f.startswith(".")
            and os.path.isfile(os.path.join(table_folder, f))
        )
    ]
    if not file_list:
        _LOGGER.warning(
            f"No CSV files found in '{table_folder}'. "
            f"Make sure Part 2 was completed successfully."
        )
        return

    summary_data = []

    for file_name in file_list:
        key = file_name.removesuffix("_orientation_distribution.csv")
        name = z_stacks_info.get(key, {}).get("original_filename", file_name)
        run = active_run()
        if run is not None and not run.eligible(name):
            continue
        with image_attempt(name, "Summary", final=True):
            file_path = os.path.join(table_folder, file_name)
            _LOGGER.info(f"\nProcessing CSV file: {file_name}")

            # Read CSV file into DataFrame
            read_file = pd.read_csv(file_path)

            # Rename columns
            read_file.rename(
                columns={
                    read_file.columns[0]: "ori_angle",
                    read_file.columns[1]: "occ_value",
                },
                inplace=True,
            )

            # Find the angle of maximum occupancy value
            angle_of_max_occ_value = read_file["ori_angle"][
                read_file["occ_value"].idxmax()
            ]

            # Normalize angles relative to the angle of maximum value
            read_file["angles_normalized_to_angle_of_MOV"] = (
                read_file["ori_angle"] - angle_of_max_occ_value
            )
            read_file["corrected_angles"] = read_file[
                "angles_normalized_to_angle_of_MOV"
            ].apply(correct_angle)

            # Rank corrected angles
            read_file["rank_of_angle_occ_value"] = read_file[
                "corrected_angles"
            ].rank(method="min")

            # Compute percentages
            sum_of_occ_values = read_file["occ_value"].sum()
            read_file["perc_occvalue2sum_of_occvalue"] = (
                read_file["occ_value"] / sum_of_occ_values
            ) * 100

            # Filter rows by angle range
            filtered_data = read_file[
                (read_file["corrected_angles"] >= -angle_value)
                & (read_file["corrected_angles"] <= angle_value)
            ]
            percentage_of_fibers_aligned_within_angle = filtered_data[
                "perc_occvalue2sum_of_occvalue"
            ].sum()

            # Determine orientation mode
            orientation_mode = "disorganized"
            if percentage_of_fibers_aligned_within_angle >= 55:
                orientation_mode = "aligned"

            # Sort DataFrame
            read_file_sorted = read_file.sort_values(
                by="rank_of_angle_occ_value"
            )

            # Generate output file name
            output_file_name = (
                f"{os.path.splitext(file_name)[0]}_processed.csv"
            )
            output_file_path = os.path.join(analysis_folder, output_file_name)
            read_file_sorted.to_csv(output_file_path, index=False)
            _LOGGER.info(f"Processed data saved at: {output_file_path}")

            processed_base_name = os.path.splitext(
                file_name.replace("_orientation_distribution.csv", "")
            )[0]

            # Get Z-stack info
            z_stacks_info_entry = z_stacks_info.get(processed_base_name, None)
            if z_stacks_info_entry is not None:
                number_of_z_stacks = z_stacks_info_entry["number_of_z_stacks"]
                z_stack_type = z_stacks_info_entry["z_stack_type"]
            else:
                number_of_z_stacks = "N/A"
                z_stack_type = "N/A"

            summary_data.append(
                {
                    "File_Name": file_name,
                    "Number_of_Z_Stacks": number_of_z_stacks,
                    "Z_Stack_Type": z_stack_type,
                    f"Percentage_Fibers_Aligned_Within_{angle_value}_Degree": (
                        percentage_of_fibers_aligned_within_angle
                    ),
                    "Orientation_Mode": orientation_mode,
                }
            )

    if not summary_data:
        return

    # Save summary data
    summary_df = pd.DataFrame(summary_data)
    summary_file_path = os.path.join(analysis_folder, "Alignment_Summary.csv")
    summary_df.to_csv(summary_file_path, index=False)
    _LOGGER.info(f"\nSummary data saved at: {summary_file_path}")

    _LOGGER.info(
        f"\nProcessing completed for folder {results_folder}. "
        f"All results saved in folder: {results_folder}"
    )


@folder_logged("folder_path")
def process_folder(
    folder_path,
    fibronectin_channel_index,
    angle_value,
    desired_width,
    desired_height,
    ij,
):
    """
    Process a single folder (Part 1, Part 2, and Part 3).

    Args:
        folder_path (str): Path to the folder for processing.
        fibronectin_channel_index (int): Index of the fibronectin
        channel.
        angle_value (float): Angle (in degrees) for analysis.
        desired_width (int): Desired width of output images.
        desired_height (int): Desired height of output images.
        ij: ImageJ context.
    """
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Source folder not found: {folder_path}")

    # Convert angle value to string for folder naming
    angle_str = f"{angle_value}".replace(".", "_")

    # Create folders for results
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S_%f")
    results_folder, table_folder, images_folder = create_results_folders(
        folder_path, angle_str, timestamp
    )

    names = image_names(folder_path, (".nd2", ".tif", ".tiff"))
    run = ImageRun(
        results_folder,
        folder_path,
        names,
        {
            "channel_index": fibronectin_channel_index,
            "angle": angle_value,
            "width": desired_width,
            "height": desired_height,
        },
        stages=PROGRESS_STAGES,
    )
    with scoped_file_log(_LOGGER, Path(results_folder)), run:
        if not names:
            return run.finish()
        run.reject_collisions()
        _LOGGER.info("Results will be saved in: %s", results_folder)
        try:
            # Part 1: Create 2D projections
            z_stacks_info_folder = process_part1(
                folder_path,
                results_folder,
                fibronectin_channel_index,
                desired_width,
                desired_height,
                ij,
            )

            # Part 2: Orientation analysis
            process_part2_orientationpy(
                results_folder, images_folder, z_stacks_info_folder
            )

            # Part 3: Process CSV files and generate summary
            analysis_folder = os.path.join(results_folder, "Analysis")
            if not os.path.exists(analysis_folder):
                os.makedirs(analysis_folder)

            process_part3(
                results_folder,
                analysis_folder,
                angle_value,
                z_stacks_info_folder,
            )
            return run.finish(Path(analysis_folder) / "Alignment_Summary.csv")
        except CANCELLATIONS:
            outcome("CANCELLED", "Alignment interrupted")
            raise
        except Exception:
            _LOGGER.exception("Alignment folder analysis failed.")
            raise


def main_fibronectin_processing(
    input_file_path, angle_value=15, desired_width=500, desired_height=500
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
    folder_paths = get_folder_paths(input_file_path)
    register_sources(folder_paths, input_file_path)
    fibr_chan_index = ask_channel()
    confirm_start()
    failed = 0
    ij = None
    try:
        for folder_path in folder_paths:
            try:
                with folder_scope(folder_path):
                    if image_names(folder_path, (".nd2", ".tif", ".tiff")):
                        if ij is None:
                            phase("Starting Fiji")
                            ij = initialize_imagej()
                    state = process_folder(
                        folder_path,
                        fibr_chan_index,
                        angle_value,
                        desired_width,
                        desired_height,
                        ij,
                    )
                    failed += state != "SUCCESS"
            except CANCELLATIONS:
                raise
            except Exception as error:
                folder_error(folder_path, error)
                failed += 1
    finally:
        if ij is not None:
            phase("Closing ImageJ context")
            ij.dispose()
    outcome(
        "FINISHED",
        f"Alignment: {len(folder_paths) - failed} complete; "
        f"{failed} incomplete or failed folder(s).",
    )
    return int(bool(failed))
