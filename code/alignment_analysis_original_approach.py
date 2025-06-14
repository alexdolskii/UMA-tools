#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script performs fibronectin orientation analysis using ImageJ and OrientationJ.
It processes microscopy images through three main steps:
1) Creates 2D projections of fibronectin channels
2) Applies orientation analysis using OrientationJ plugin
3) Processes results and generates alignment statistics
"""

import os
import sys
import json
import time
from datetime import datetime
from pathlib import Path

import imagej
import pandas as pd
from scyjava import jimport


class ImageJInitializationError(Exception):
    """Exception raised for errors during ImageJ initialization."""
    pass


class JSONValidationError(Exception):
    """Exception raised for invalid JSON structure or content."""
    pass


def get_folder_paths(json_path):
    """
    Read folder paths from JSON file and validate structure.
    
    Args:
        json_path (str): Path to JSON file
        
    Returns:
        list: Valid folder paths
        
    Raises:
        JSONValidationError: If JSON structure is invalid
        FileNotFoundError: If JSON file doesn't exist
    """
    print(f"\nLooking for input_paths.json at: {json_path}")
    
    if not os.path.isfile(json_path):
        raise FileNotFoundError(f"input_paths.json not found at {json_path}")

    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_data = json.load(f)
        
        # Validate JSON structure
        if "folder_paths" not in json_data:
            raise JSONValidationError("JSON missing 'folder_paths' key")
            
        folder_paths = json_data["folder_paths"]
        
        if not isinstance(folder_paths, list):
            raise JSONValidationError("'folder_paths' should be a list")
            
        if not folder_paths:
            raise JSONValidationError("'folder_paths' list is empty")
            
        return folder_paths
        
    except Exception as e:
        raise JSONValidationError(f"Error reading input_paths.json: {e}")


def get_fibronectin_channel(folder_paths):
    """
    Prompt user for fibronectin channel configuration.
    
    Args:
        folder_paths (list): List of folder paths to process
        
    Returns:
        dict: Fibronectin channel indices for each folder
    """
    fibronectin_channel_indices = {}
    
    same_channel = input("\nDo all folders have the same fibronectin channel? (yes/no): ").strip().lower()
    
    if same_channel in ('yes', 'y'):
        channel_input = input("Enter the fibronectin channel number (starting from 1): ").strip()
        try:
            fibronectin_channel_index = int(channel_input) - 1
            if fibronectin_channel_index < 0:
                print("Using default channel 1")
                fibronectin_channel_index = 0
        except ValueError:
            print("Invalid input. Using default channel 1.")
            fibronectin_channel_index = 0
        
        for folder in folder_paths:
            fibronectin_channel_indices[folder] = fibronectin_channel_index
    else:
        for folder in folder_paths:
            print(f"\nFolder: {folder}")
            channel_input = input("Enter the fibronectin channel number (starting from 1): ").strip()
            try:
                fibronectin_channel_index = int(channel_input) - 1
                if fibronectin_channel_index < 0:
                    print("Using default channel 1")
                    fibronectin_channel_index = 0
            except ValueError:
                print("Invalid input. Using default channel 1.")
                fibronectin_channel_index = 0
            fibronectin_channel_indices[folder] = fibronectin_channel_index
    
    return fibronectin_channel_indices


def create_results_folders(input_folder, angle_str, timestamp):
    """
    Create folder structure for analysis results.
    
    Args:
        input_folder (str): Path to input folder
        angle_str (str): Formatted angle value for folder naming
        timestamp (str): Timestamp for folder naming
        
    Returns:
        tuple: Paths to (results_folder, excel_folder, images_folder)
    """
    results_folder_name = f"Alignment_assay_results_angle_{angle_str}_{timestamp}"
    results_folder = os.path.join(input_folder, results_folder_name)
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    
    excel_folder = os.path.join(results_folder, "Excel")
    images_folder = os.path.join(results_folder, "Images")
    Path(excel_folder).mkdir(parents=True, exist_ok=True)
    Path(images_folder).mkdir(parents=True, exist_ok=True)
    
    return results_folder, excel_folder, images_folder


def process_image_file(ij, file_path, fibronectin_channel_index, desired_width, desired_height):
    """
    Process an individual image file through projection and channel extraction.
    
    Args:
        ij: ImageJ instance
        file_path (str): Path to image file
        fibronectin_channel_index (int): Channel index to extract
        desired_width (int): Target output width
        desired_height (int): Target output height
        
    Returns:
        tuple: (processed_imp, z_stack_info)
    """
    IJ = jimport('ij.IJ')
    ZProjector = jimport('ij.plugin.ZProjector')
    ImagePlus = jimport('ij.ImagePlus')
    
    try:
        # Open image using Bio-Formats
        img = ij.io().open(file_path)
        if img is None:
            print(f"Failed to open image: {file_path}")
            return None, None
            
        # Convert to ImagePlus
        imp = ij.convert().convert(img, ImagePlus)
        dimensions = imp.getDimensions()
        width, height, channels, slices, frames = dimensions
        
        # Determine initial Z-stack info
        initial_z_stacks = slices if slices > 1 else channels
        z_stack_type = 'slices' if slices > 1 else 'channels'
        
        # Process based on file type
        file_ext = os.path.splitext(file_path)[1].lower()
        
        if file_ext in ('.nd2', '.oif', '.oib'):
            # Extract fibronectin channel
            if channels > 1:
                imp.setC(fibronectin_channel_index + 1)
                IJ.run(imp, "Make Substack...", f"channels={fibronectin_channel_index+1} slices=1-{slices} frames=1-{frames}")
                imp.close()
                imp = IJ.getImage()
                channels = imp.getNChannels()
        
        elif file_ext in ('.tif', '.tiff'):
            # Project across channels if needed
            if channels > 1:
                imp.show()
                IJ.run("Z Project...", "projection=[Max Intensity] all")
                imp.close()
                imp = IJ.getImage()
                channels = imp.getNChannels()
        
        # Project along slices if needed
        slices = imp.getNSlices()
        if slices > 1:
            zp = ZProjector(imp)
            zp.setMethod(ZProjector.MAX_METHOD)
            zp.doProjection()
            projected_imp = zp.getProjection()
            imp.close()
            imp = projected_imp
        
        # Resize image
        IJ.run(imp, "Size...", f"width={desired_width} height={desired_height} constrain average interpolation=Bilinear")
        
        return imp, {
            'number_of_z_stacks': initial_z_stacks,
            'z_stack_type': z_stack_type
        }
        
    except Exception as e:
        print(f"Error processing image: {e}")
        return None, None


def process_part1(ij, folder_path, fibronectin_channel_index, results_folder, desired_width, desired_height):
    """
    Part 1: Create 2D projections of fibronectin channels.
    
    Args:
        ij: ImageJ instance
        folder_path (str): Input folder path
        fibronectin_channel_index (int): Channel index to extract
        results_folder (str): Output folder for results
        desired_width (int): Target image width
        desired_height (int): Target image height
        
    Returns:
        dict: Z-stack information for processed files
    """
    print(f"\nProcessing folder: {folder_path}")
    IJ = jimport('ij.IJ')
    z_stacks_info = {}
    
    for filename in os.listdir(folder_path):
        if filename.startswith('._'):
            continue
            
        if not filename.lower().endswith(('.tif', '.tiff', '.nd2', '.oif', '.oib')):
            continue
            
        file_path = os.path.join(folder_path, filename)
        print(f"\nProcessing file: {file_path}")
        
        imp, z_stack_info = process_image_file(
            ij, 
            file_path, 
            fibronectin_channel_index, 
            desired_width, 
            desired_height
        )
        
        if imp is None:
            continue
            
        # Save processed image
        output_filename = os.path.splitext(filename)[0] + '_processed.tif'
        output_path = os.path.join(results_folder, output_filename)
        IJ.saveAs(imp, "Tiff", output_path)
        imp.close()
        
        # Store Z-stack info using ORIGINAL base name as key
        original_base_name = os.path.splitext(filename)[0]
        z_stacks_info[original_base_name] = {
            'original_filename': filename,
            'number_of_z_stacks': z_stack_info['number_of_z_stacks'],
            'z_stack_type': z_stack_info['z_stack_type']
        }
        
    print(f"Part 1 completed for folder {folder_path}")
    return z_stacks_info

def process_part2(ij, results_folder, images_folder):
    """
    Part 2: Apply OrientationJ analysis to processed images.
    
    Args:
        ij: ImageJ instance
        results_folder (str): Folder with processed images
        images_folder (str): Output folder for analysis images
    """
    print("\nStarting Part 2: Applying OrientationJ Plugin...")
    IJ = jimport('ij.IJ')
    WindowManager = jimport('ij.WindowManager')
    
    processed_files = [f for f in os.listdir(results_folder) 
                      if f.endswith('_processed.tif') and not f.startswith('._')]
    
    for filename in processed_files:
        file_path = os.path.join(results_folder, filename)
        print(f"\nProcessing file: {file_path}")
        
        try:
            # Open image
            imp = IJ.openImage(file_path)
            if imp is None:
                continue
                
            # OrientationJ Analysis
            imp.show()
            IJ.run("OrientationJ Analysis", "tensor=3.0 gradient=4 color-survey=on hsb=on hue=Orientation sat=Coherency bri=Original-Image radian=on")
            time.sleep(0.5)
            
            # Save analysis image
            analysis_title = "OJ-Color-survey-1"
            analysis_imp = WindowManager.getImage(analysis_title)
            if analysis_imp:
                analysis_filename = filename.replace('_processed.tif', '_oj_analysis.tif')
                analysis_path = os.path.join(images_folder, analysis_filename)
                IJ.saveAs(analysis_imp, "Tiff", analysis_path)
                analysis_imp.close()
                
            # Close all and reopen for distribution
            IJ.run("Close All")
            imp = IJ.openImage(file_path)
            imp.show()
            
            # OrientationJ Distribution
            IJ.run("OrientationJ Distribution", "tensor=3.0 gradient=4 radian=on histogram=on table=on min-coherency=0.0 min-energy=0.0")
            time.sleep(0.5)
            
            # Save results table
            excel_filename = filename.replace('_processed.tif', '_oj_distribution.csv')
            excel_path = os.path.join(os.path.dirname(images_folder), "Excel", excel_filename)
            IJ.saveAs("Results", excel_path)
            IJ.run("Clear Results")
            
            # Clean up
            IJ.run("Close All")
            
        except Exception as e:
            print(f"Error during OrientationJ processing: {e}")
            IJ.run("Close All")


def process_csv_file(file_path, angle_value):
    """
    Process an OrientationJ CSV result file.
    
    Args:
        file_path (str): Path to CSV file
        angle_value (float): Alignment angle threshold
        
    Returns:
        tuple: (processed_df, alignment_percentage, orientation_mode)
    """
    df = pd.read_csv(file_path)
    
    # Rename columns
    df.columns = ['orientation_angle', 'occurrence_value']
    
    # Find peak angle
    max_occurrence_idx = df['occurrence_value'].idxmax()
    angle_of_max = df.loc[max_occurrence_idx, 'orientation_angle']
    
    # Normalize angles
    df['angles_normalized_to_angle_of_MOV'] = df['orientation_angle'] - angle_of_max
    df['corrected_angles'] = df['angles_normalized_to_angle_of_MOV'].apply(
        lambda x: x + 180 if x < -90 else (x - 180 if x > 90 else x)
    )
    
    # Calculate statistics
    total_occurrence = df['occurrence_value'].sum()
    df['percent_occurrence'] = (df['occurrence_value'] / total_occurrence) * 100
    
    # Calculate alignment percentage
    aligned_mask = (df['corrected_angles'] >= -angle_value) & (df['corrected_angles'] <= angle_value)
    alignment_percentage = df.loc[aligned_mask, 'percent_occurrence'].sum()
    
    # Determine orientation mode
    orientation_mode = 'aligned' if alignment_percentage >= 55 else 'disorganized'
    
    return df, alignment_percentage, orientation_mode


def process_part3(results_folder, analysis_folder, angle_value, z_stacks_info, timestamp):
    """
    Part 3: Process OrientationJ results and generate summary.
    
    Args:
        results_folder (str): Folder with analysis results
        analysis_folder (str): Output folder for processed data
        angle_value (float): Alignment angle threshold
        z_stacks_info (dict): Z-stack information from Part 1
        timestamp (str): Timestamp for file naming
    """
    print("\nStarting Part 3: Processing Excel Files...")
    excel_folder = os.path.join(results_folder, "Excel")
    
    if not os.path.exists(excel_folder):
        print("Excel folder not found. Skipping Part 3.")
        return
        
    csv_files = [f for f in os.listdir(excel_folder) 
                if f.endswith('.csv') and not f.startswith('._')]
    
    summary_data = []
    
    for csv_file in csv_files:
        csv_path = os.path.join(excel_folder, csv_file)
        print(f"\nProcessing: {csv_file}")
        
        try:
            # Process CSV file
            df, alignment_percentage, orientation_mode = process_csv_file(csv_path, angle_value)
            
            # Save processed data
            output_filename = csv_file.replace('.csv', f'_processed_{timestamp}.xlsx')
            output_path = os.path.join(analysis_folder, output_filename)
            df.to_excel(output_path, index=False)
            
            # Get ORIGINAL base name from CSV filename
            # Example: 'experiment1_oj_distribution.csv' -> 'experiment1'
            base_name = csv_file.replace('_oj_distribution.csv', '')
            base_name = base_name.replace('_processed', '')  # Remove processed suffix
            
            # Get Z-stack info using original base name
            z_info = z_stacks_info.get(base_name, {})
            
            # Add to summary
            summary_data.append({
                'File_Name': csv_file,
                'Number_of_Z_Stacks': z_info.get('number_of_z_stacks', 'N/A'),
                'Z_Stack_Type': z_info.get('z_stack_type', 'N/A'),
                f'Percentage_Fibers_Aligned_Within_{angle_value}_Degree': alignment_percentage,
                'Orientation_Mode': orientation_mode
            })
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
    
    # Save summary
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(analysis_folder, f'Alignment_Summary_{timestamp}.xlsx')
        summary_df.to_excel(summary_path, index=False)
        print(f"\nSummary saved to: {summary_path}")

def process_folder(
    ij, 
    folder_path, 
    fibronectin_channel_index, 
    angle_value, 
    desired_width, 
    desired_height
):
    """
    Process a single folder through all analysis steps.
    
    Args:
        ij: ImageJ instance
        folder_path (str): Path to folder to process
        fibronectin_channel_index (int): Fibronectin channel index
        angle_value (float): Alignment angle threshold
        desired_width (int): Target image width
        desired_height (int): Target image height
    """
    # Create results folders
    angle_str = f"{angle_value}".replace('.', '_')
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_folder, excel_folder, images_folder = create_results_folders(
        folder_path, angle_str, timestamp
    )
    
    # Part 1: Create 2D projections
    z_stacks_info = process_part1(
        ij, 
        folder_path, 
        fibronectin_channel_index, 
        results_folder, 
        desired_width, 
        desired_height
    )
    
    # Part 2: OrientationJ analysis
    process_part2(ij, results_folder, images_folder)
    
    # Part 3: Process results
    analysis_folder = os.path.join(results_folder, 'Analysis')
    os.makedirs(analysis_folder, exist_ok=True)
    
    process_part3(
        results_folder,
        analysis_folder,
        angle_value,
        z_stacks_info,
        timestamp
    )
    
    normalize_saved_images(results_folder)

    print(f"\nProcessing completed for folder: {folder_path}")
    print(f"Results saved in: {results_folder}")


def main():
    """Main function to execute the fibronectin analysis workflow."""
    try:
        # Initialize ImageJ
        script_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(script_dir)
        fiji_path = os.path.join(parent_dir, "Fiji.app")
        
        if not os.path.exists(fiji_path):
            raise FileNotFoundError(f"Fiji.app not found at {fiji_path}")
        
        print(f"Initializing ImageJ from: {fiji_path}")
        ij = imagej.init(fiji_path, mode='interactive')
        print("ImageJ initialization completed.")
        
        # Get folder paths
        json_path = os.path.join(parent_dir, "input_paths.json")
        folder_paths = get_folder_paths(json_path)
        print(f"Found {len(folder_paths)} folders to process")
        
        # Get angle value
        angle_input = input("\nEnter analysis angle (default 15): ").strip()
        angle_value = float(angle_input) if angle_input else 15.0
        angle_str = f"{angle_value}".replace('.', '_')
        
        # Get fibronectin channels
        fibronectin_channel_indices = get_fibronectin_channel(folder_paths)
        
        # Confirm processing
        start_processing = input("\nStart processing? (yes/no): ").strip().lower()
        if start_processing not in ('yes', 'y'):
            print("Processing canceled by user.")
            return
            
        # Process each folder
        desired_width, desired_height = 500, 500
        for folder_path in folder_paths:
            process_folder(
                ij,
                folder_path,
                fibronectin_channel_indices[folder_path],
                angle_value,
                desired_width,
                desired_height
            )
            
        print("\nAll folders processed successfully.")
        
    except Exception as e:
        print(f"\nError: {e}")
        sys.exit(1)
        
    finally:
        # Clean up ImageJ
        if 'ij' in locals():
            print("Disposing ImageJ context...")
            ij.dispose()
        print("Script execution completed.")

import numpy as np
import matplotlib.colors
from skimage import io, img_as_ubyte
import matplotlib.pyplot as plt

def normalize_saved_images(results_folder):
    """
    Part 4: Normalize saved orientation images using angle data from CSV files.
    
    Args:
        results_folder (str): Main results folder containing Excel and Images subfolders
    """
    print("\nStarting Part 4: Normalizing saved orientation images...")
    
    # Create normalized images subfolder
    images_folder = os.path.join(results_folder, "Images")
    normalized_folder = os.path.join(images_folder, "normalized_images")
    Path(normalized_folder).mkdir(parents=True, exist_ok=True)
    
    excel_folder = os.path.join(results_folder, "Excel")
    
    # Get all distribution CSV files
    csv_files = [f for f in os.listdir(excel_folder) 
                if f.endswith('_oj_distribution.csv') and not f.startswith('._')]
    
    if not csv_files:
        print(f"No distribution CSV files found in {excel_folder}")
        return
    
    for csv_file in csv_files:
        try:
            print(f"\nProcessing {csv_file}")
            
            # Extract base filename from CSV name
            base_name = csv_file.replace('_oj_distribution.csv', '')
            image_filename = base_name + '_oj_analysis.tif'
            image_path = os.path.join(images_folder, image_filename)
            
            print(f"Looking for image: {image_path}")
            
            if not os.path.exists(image_path):
                print(f"Image not found: {image_path}")
                continue
                
            # Read CSV file
            csv_path = os.path.join(excel_folder, csv_file)
            df = pd.read_csv(csv_path)
            
            # Find modal angle (max occurrence)
            max_occurrence_idx = df.iloc[:, 1].idxmax()
            modal_angle_deg = df.iloc[max_occurrence_idx, 0]
            print(f"Modal angle: {modal_angle_deg:.2f}°")
            
            # Load the orientation image
            img = io.imread(image_path)
            
            # If image is RGBA, convert to RGB
            if img.shape[-1] == 4:
                img = img[..., :3]
            
            # Convert to float and normalize to [0, 1]
            img_float = img.astype(np.float32) / 255.0
            
            # Create HSV representation
            hsv_img = matplotlib.colors.rgb_to_hsv(img_float)
            
            # Calculate hue shift: -2 * modal_angle_deg
            hue_shift = -2 * modal_angle_deg
            print(f"Applying hue shift: {hue_shift:.2f}°")
            
            # Apply hue shift (convert to [0, 1] range)
            hsv_img[..., 0] = (hsv_img[..., 0] + (hue_shift / 360)) % 1.0
            
            # Convert back to RGB
            normalized_rgb = matplotlib.colors.hsv_to_rgb(hsv_img)
            
            # Create figure for saving
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.imshow(normalized_rgb)
            ax.axis('off')
            
            # Add colorbar
            sm = plt.cm.ScalarMappable(
                cmap='hsv',
                norm=plt.Normalize(vmin=-90, vmax=90))
            sm.set_array([])
            
            cbar = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.7)
            cbar.set_label("Deviation from Dominant Direction (°)")
            
            # Save normalized image
            normalized_filename = base_name + '_oj_analysis_normalized.png'
            normalized_path = os.path.join(normalized_folder, normalized_filename)
            plt.savefig(normalized_path, bbox_inches='tight', pad_inches=0)
            plt.close(fig)
            
            print(f"Saved normalized image: {normalized_path}")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    main()