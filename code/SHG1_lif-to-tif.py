#!/usr/bin/env python3

import os
import json
import tifffile
import numpy as np
import argparse
import psutil  # For monitoring memory usage
import scyjava as sj
from datetime import datetime
from readlif.reader import LifFile

# Increase memory for JVM
os.environ['_JAVA_OPTIONS'] = (
    "-Xmx16g "  # Up to 16GB of RAM
    "-XX:+IgnoreUnrecognizedVMOptions "
    "--illegal-access=warn "
    "--add-opens=java.base/java.lang=ALL-UNNAMED "
)

import imagej  # Import ImageJ after setting the environment


class LifFileProcessorError(Exception):
    """Custom exception for LIF file processing errors."""
    pass


class ImageJInitializationError(Exception):
    """Exception for ImageJ initialization errors."""
    pass


def validate_json_input(json_path):
    """
    Reads and validates a JSON file, which must contain the key 'paths_to_files'
    with a list of folder paths.
    """
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            config = json.load(f)
        if "paths_to_files" not in config or not config["paths_to_files"]:
            raise LifFileProcessorError(
                "JSON file does not contain 'paths_to_files' key or it is empty."
            )
        return config["paths_to_files"]
    except Exception as e:
        raise LifFileProcessorError(
            f"Error reading or parsing JSON file '{json_path}': {e}"
        )


def create_output_folder(base_folder):
    """
    Creates the main output folder with the name:
      output_<folder_name>_<timestamp>
    Returns the path to this folder.
    """
    folder_name = os.path.basename(os.path.normpath(base_folder))
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_base_dir = f"output_{folder_name}_{timestamp}"
    os.makedirs(output_base_dir, exist_ok=True)
    return output_base_dir


def create_processed_folder(base_folder):
    """
    Creates a folder for processed 2D projections:
      output_processed_<folder_name>_<timestamp>
    Returns the path to this folder.
    """
    folder_name = os.path.basename(os.path.normpath(base_folder))
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    processed_dir = f"output_processed_{folder_name}_{timestamp}"
    os.makedirs(processed_dir, exist_ok=True)
    return processed_dir


def ensure_channel_dirs(output_base_dir, max_channels):
    """
    Creates subfolders channel_1, channel_2, ..., channel_N
    inside output_base_dir. Returns a list of paths to these folders.
    """
    channel_dirs = []
    for ch_index in range(1, max_channels + 1):
        dir_path = os.path.join(output_base_dir, f"channel_{ch_index}")
        os.makedirs(dir_path, exist_ok=True)
        channel_dirs.append(dir_path)
    return channel_dirs


def process_lif_file(lif_path, channel_dirs):
    """
    Processes a single LIF file:
      - Opens the file
      - Iterates through each scene
      - Saves Z-stacks as TIFF files per channel and mosaic tile
      - Saves results in channel_dirs[c]
    """
    file_prefix = os.path.splitext(os.path.basename(lif_path))[0]

    try:
        lif_file = LifFile(lif_path)
        images = list(lif_file.get_iter_image())
        if not images:
            raise LifFileProcessorError(f"No scenes found in file '{lif_path}'.")
    except Exception as e:
        raise LifFileProcessorError(f"Error opening LIF file '{lif_path}': {e}")

    print(f"\n[Processing LIF] {lif_path} — {len(images)} scene(s) found.")

    for image in images:
        scene_name = image.info.get('name', 'scene_Unknown')
        x_size, y_size, z_stacks, _, mosaic_tiles = image.dims
        channels = image.channels

        for c in range(channels):
            if c >= len(channel_dirs):
                print(f"[Warning] Skipping channel {c+1}, no corresponding folder.")
                continue

            channel_folder = channel_dirs[c]

            for m_index in range(mosaic_tiles):
                stack_3d = np.zeros((z_stacks, y_size, x_size), dtype=np.uint16)
                for z in range(z_stacks):
                    frame = image.get_frame(z=z, t=0, c=c, m=m_index)
                    stack_3d[z] = frame

                safe_scene_name = scene_name.replace(" ", "_")
                out_name = f"{file_prefix}_{safe_scene_name}_tile_{m_index+1}_C{c+1}.tif"
                out_path = os.path.join(channel_folder, out_name)

                tifffile.imwrite(out_path, stack_3d, photometric="minisblack")
                print(f"  Saved stack: {out_path}")


def initialize_imagej():
    """
    Initializes ImageJ2 in headless mode.
    """
    print("\nInitializing ImageJ2 in headless mode...")
    try:
        ij = imagej.init('sc.fiji:fiji', mode='headless')
    except Exception as e:
        raise ImageJInitializationError(f"Error initializing ImageJ: {e}")
    print("ImageJ2 successfully initialized.")
    return ij


def main():
    parser = argparse.ArgumentParser(
        description="LIF file processing and Z-projection using ImageJ."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Path to JSON file containing 'paths_to_files' key.")
    args = parser.parse_args()

    try:
        folder_paths = validate_json_input(args.input)
    except LifFileProcessorError as e:
        print(f"[ERROR] {e}")
        return

    try:
        ij = initialize_imagej()
    except ImageJInitializationError as e:
        print(f"[ERROR] {e}")
        return

    for folder in folder_paths:
        if not os.path.isdir(folder):
            print(f"[WARNING] Folder '{folder}' does not exist, skipping.")
            continue

        print(f"\n=== Processing folder: {folder} ===")

        output_base_dir = create_output_folder(folder)
        output_processed_dir = create_processed_folder(folder)

        lif_files = [f for f in os.listdir(folder) if f.lower().endswith(".lif")]
        print(f"  Found {len(lif_files)} LIF file(s).")

        if not lif_files:
            print("  No LIF files, skipping.")
            continue

        global_max_channels = max(
            (img.channels for lf in lif_files for img in LifFile(os.path.join(folder, lf)).get_iter_image()),
            default=0
        )

        if global_max_channels == 0:
            print("  No valid channels, skipping.")
            continue

        ensure_channel_dirs(output_base_dir, global_max_channels)
        print("\nFolder processing completed.")

    print("\nAll folders processed.")


if __name__ == "__main__":
    main()
