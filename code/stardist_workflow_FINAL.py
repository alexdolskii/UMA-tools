#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Nuclei Analysis Pipeline
- Processes ND2 files to extract the nuclei channel
- Resizes images for StarDist 3D
- Performs StarDist 3D segmentation
- Analyzes masks for nuclear features
- Performs clustering analysis to identify nuclear layers
- Outputs CSV with quantitative results and QC images
"""

import os
import json
import argparse
import numpy as np
import matplotlib
matplotlib.rcParams["image.interpolation"] = 'none'
import matplotlib.pyplot as plt
import pandas as pd
import logging
from datetime import datetime
from pathlib import Path
from tifffile import imread, imwrite
from skimage.color import label2rgb
from skimage.measure import regionprops_table
from skimage.filters import gaussian
from skimage.filters import median
from scipy.ndimage import uniform_filter
from skimage.morphology import ball
from csbdeep.utils import normalize
from stardist.models import StarDist3D
import imagej
from scyjava import jimport
from sklearn.cluster import DBSCAN, HDBSCAN
from sklearn.datasets import make_blobs
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
from collections import Counter

# Initialize logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('nuclei_analysis.log'),
        logging.StreamHandler()
    ]
)

np.random.seed(6)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='3D Nuclei Analysis Pipeline')
    parser.add_argument('-i', '--input', type=str, required=False,
                        help='Path to the input configuration JSON file',
                        default="input_paths.json")
    return parser.parse_args()

def load_config(config_path):
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Configuration file '{config_path}' not found.")
    with open(config_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Extract folder paths
    folders = data.get('folder_paths', [])
    valid_folders = [f for f in folders if os.path.isdir(f)]
    if not valid_folders:
        raise ValueError("No valid folders provided in configuration.")
    
    # Extract other configuration parameters
    config = {
        'folders': valid_folders,
        'nuclei_channel': data.get('nuclei_channel', 1),  # Default to channel 1 if not specified
        'gaussian_sigma': data.get('gaussian_sigma', 4.0),
        'mean_radius': data.get('mean_radius', 3),
        'z_scale_factor': data.get('z_scale_factor', 32),
        'model_path': data.get('model_path')  # extract model_path
    }
    
    return config

def extract_and_resize_nuclei(input_dir, output_dir, nuclei_channel=1, target_size=(1024, 1024), gaussian_sigma=4.0, mean_radius=3):
    """
    Extract and resize the nuclei channel from ND2 files.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing ND2 files
    output_dir : str
        Directory to save processed TIF files
    nuclei_channel : int
        1-based channel number for nuclei (default=1)
    target_size : tuple
        Target size for resizing (width, height)
    gaussian_sigma : float
        Sigma value for Gaussian filter
    mean_radius : int
        Radius for Mean filter
    """
    # Initialize ImageJ (headless)
    ij = imagej.init('sc.fiji:fiji', mode='headless')
    IJ = jimport('ij.IJ')
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"Using nuclei channel = {nuclei_channel}")

    for file in os.listdir(input_dir):
        if not file.lower().endswith('.nd2') or file.startswith(('.', '_')):
            continue

        file_path = os.path.join(input_dir, file)
        output_path = os.path.join(output_dir, f"{Path(file).stem}_nuclei_resized.tif")
        print(f"Processing: {file}")

        imp = IJ.openImage(file_path)
        if imp is None:
            logging.warning(f"Could not open: {file_path}")
            continue

        width, height, channels, slices, frames = imp.getDimensions()
        if not (1 <= nuclei_channel <= channels):
            logging.error(f"Invalid channel {nuclei_channel} for {file}. File has {channels} channels.")
            imp.close()
            continue

        # Extract just the nuclei channel
        IJ.run(imp, "Duplicate...", f"title=imp_nuclei duplicate channels={nuclei_channel}")
        imp_nuclei = IJ.getImage()
        
        # Convert to 8-bit before applying filters
        print("Converting to 8-bit...")
        IJ.run(imp_nuclei, "8-bit", "")
        
        # Resize the image
        IJ.run(imp_nuclei, "Size...", f"width={target_size[1]} height={target_size[0]} depth={slices} interpolation=Bilinear")
        
        # Apply Gaussian filter directly in ImageJ
        print("Applying Gaussian filter...")
        IJ.run(imp_nuclei, "Gaussian Blur 3D...", f"x={gaussian_sigma} y={gaussian_sigma} z={gaussian_sigma}")
        
        # Apply Mean filter directly in ImageJ
        print("Applying Mean filter...")
        IJ.run(imp_nuclei, "Mean 3D...", f"x={mean_radius} y={mean_radius} z={mean_radius}")

        # Apply Gaussian filter directly in ImageJ
        print("Applying Gaussian filter...")
        IJ.run(imp_nuclei, "Gaussian Blur 3D...", f"x={gaussian_sigma} y={gaussian_sigma} z={gaussian_sigma}")
        
        # Check the pixel values after filtering
        stack = imp_nuclei.getStack()
        middle_slice = slices // 2
        processor = stack.getProcessor(middle_slice + 1)  # 1-based indexing
        sample_pixels = processor.getPixels()
        try:
            # Convert Java array to Python list and print a sample
            sample_list = [sample_pixels[i] for i in range(min(10, len(sample_pixels)))]
            print(f"Sample pixels from middle slice after filtering: {sample_list}")
            
            # Get min and max values
            stats = processor.getStatistics()
            print(f"Min: {stats.min}, Max: {stats.max}, Mean: {stats.mean}")
        except Exception as e:
            print(f"Error getting pixel statistics: {e}")
        
        # Save the filtered image
        IJ.saveAs(imp_nuclei, "Tiff", output_path)
        imp_nuclei.close()
        imp.close()
        print(f"Saved filtered image: {output_path}")

def run_stardist_segmentation(folders, model):
    for folder in folders:
        processed_dir = os.path.join(folder, "processed")
        masks_dir = os.path.join(folder, "masks")
        os.makedirs(masks_dir, exist_ok=True)

        tiff_files = [f for f in os.listdir(processed_dir) 
              if f.lower().endswith(('.tif', '.tiff')) 
              and not f.startswith(('.', '_'))]
        if not tiff_files:
            print(f"No TIFF files in {processed_dir}.")
            continue

        print(f"Running StarDist on {len(tiff_files)} files in {processed_dir}...")

        for file in tiff_files:
            try:
                # Load image and ensure it's in the correct format for StarDist
                img_path = os.path.join(processed_dir, file)
                img = imread(img_path)
                
                # Normalize the image
                img_norm = normalize(img, 1, 99.8, axis=None)
                
                # Make sure the image has the correct dimensions
                if img_norm.ndim == 2:
                    img_norm = img_norm[np.newaxis, ...]
                
                # Print image shape for debugging
                print(f"Image shape: {img_norm.shape}")
                
                # Predict using the model
                labels, details = model.predict_instances(img_norm)
                
                base_name = Path(file).stem
                imwrite(os.path.join(masks_dir, f"{base_name}_mask.tif"), labels)

                # Create a QC image using a middle Z-slice
                mid_z = img_norm.shape[0] // 2
                plt.figure(figsize=(16, 8))
                plt.subplot(121)
                plt.imshow(img_norm[mid_z], cmap='gray')
                plt.title('Original (mid-Z)')
                plt.axis('off')
                plt.subplot(122)
                plt.imshow(img_norm[mid_z], cmap='gray')
                plt.imshow(labels[mid_z], cmap='jet', alpha=0.5)
                plt.title('Mask Overlay')
                plt.axis('off')
                plt.savefig(os.path.join(masks_dir, f"{base_name}_QC.png"), dpi=150, bbox_inches='tight')
                plt.close()

                print(f"Processed {file}: {len(np.unique(labels)) - 1} nuclei.")
            except Exception as e:
                print(f"Error with {file}: {e}")
                logging.exception(f"Segmentation failed for {file}: {str(e)}")

def analyze_masks(mask_path, output_csv, min_volume=0.5, original_size=(1024, 1024)):
    masks = imread(mask_path)
    print(f"Analyzing: {mask_path}, shape: {masks.shape}")

    scale_y = original_size[0] / masks.shape[1]
    scale_x = original_size[1] / masks.shape[2]

    props = regionprops_table(
        masks,
        properties=(
            'label', 'area', 'centroid', 'bbox', 'solidity', 'equivalent_diameter_area'
        )
    )
    df = pd.DataFrame(props)
    df['area'] *= scale_y * scale_x
    df['centroid-1'] *= scale_y
    df['centroid-2'] *= scale_x

    voxel_volume = (0.1 * scale_x) * (0.1 * scale_y) * 0.3
    df['volume_um3'] = df['area'] * voxel_volume
    df = df[df['area'] >= min_volume]

    df = df.rename(columns={
        'label': 'nucleus_id',
        'area': 'volume_voxels',
        'centroid-0': 'z_pos',
        'centroid-1': 'y_pos',
        'centroid-2': 'x_pos',
        'equivalent_diameter_area': 'diameter_um'
    })

    df.to_csv(output_csv, index=False)
    generate_qc_plots(masks, df, output_csv)
    return df

def generate_qc_plots(mask_data, df, output_csv):
    plt.figure(figsize=(20, 12))
    max_label = mask_data.max()
    rng = np.random.RandomState(42)
    label_colors = rng.rand(max_label+1, 3)

    shape_z, shape_y, shape_x = mask_data.shape
    rgb_data = np.zeros((shape_z, shape_y, shape_x, 3), dtype=np.float32)
    for lbl in range(1, max_label+1):
        rgb_data[mask_data == lbl] = label_colors[lbl]

    coverage_mask = np.any(mask_data != 0, axis=0)
    coverage_percent = 100 * coverage_mask.sum() / coverage_mask.size

    xy = rgb_data.max(axis=0)
    xz = np.flipud(rgb_data.max(axis=1))
    yz = np.flipud(rgb_data.max(axis=2))

    # Projections
    plt.subplot(2, 3, 1)
    plt.imshow(xy)
    plt.title(f'XY Projection ({len(df)} nuclei)\nCoverage: {coverage_percent:.2f}%')
    plt.xlabel('X →')
    plt.ylabel('Y →')

    plt.subplot(2, 3, 2)
    plt.imshow(xz)
    plt.title('XZ Projection')
    plt.xlabel('X →')
    plt.ylabel('Z (Top) ↓')

    plt.subplot(2, 3, 3)
    plt.imshow(yz)
    plt.title('YZ Projection')
    plt.xlabel('Y →')
    plt.ylabel('Z (Top) ↓')

    # Centroids
    for idx, (x, y, label) in enumerate(zip(df['x_pos'], df['y_pos'], df['nucleus_id'])):
        plt.subplot(2, 3, 4)
        plt.scatter(x, y, s=df.iloc[idx]['volume_voxels']/50, c=[label_colors[label]], alpha=0.7, edgecolor='k')
    plt.title('XY Centroids')
    plt.xlabel('X →')
    plt.ylabel('Y →')
    plt.grid(True)

    for idx, (x, z, label) in enumerate(zip(df['x_pos'], df['z_pos'], df['nucleus_id'])):
        plt.subplot(2, 3, 5)
        z_flipped = shape_z - z
        plt.scatter(x, z_flipped, s=df.iloc[idx]['volume_voxels']/50, c=[label_colors[label]], alpha=0.7, edgecolor='k')
    plt.title('XZ Centroids')
    plt.xlabel('X →')
    plt.ylabel('Z (Top) ↓')
    plt.grid(True)

    for idx, (y, z, label) in enumerate(zip(df['y_pos'], df['z_pos'], df['nucleus_id'])):
        plt.subplot(2, 3, 6)
        z_flipped = shape_z - z
        plt.scatter(y, z_flipped, s=df.iloc[idx]['volume_voxels']/50, c=[label_colors[label]], alpha=0.7, edgecolor='k')
    plt.title('YZ Centroids')
    plt.xlabel('Y →')
    plt.ylabel('Z (Top) ↓')
    plt.grid(True)

    plt.tight_layout()
    plot_path = os.path.splitext(output_csv)[0] + '_3D_QC.png'
    plt.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved QC plots to: {plot_path}")

def cluster_nuclei(analysis_csv_path, output_dir, z_scale_factor=32, summary_data=None):
    """
    Perform HDBSCAN clustering on nuclei coordinates from analysis CSV files
    and save the resulting plots (XZ and XY projections).
    
    Parameters:
    -----------
    analysis_csv_path : str
        Path to the analysis CSV file with nuclei coordinates
    output_dir : str
        Directory to save output plots
    z_scale_factor : int
        Factor to scale Z coordinates by (default=32)
    summary_data : dict, optional
        Dictionary to collect summary data for multiple images
    
    Returns:
    --------
    tuple
        (cluster_results DataFrame, summary_dict with counts per layer)
    """
    print(f"Running clustering analysis on: {analysis_csv_path}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load nuclei data
    data_3d = pd.read_csv(analysis_csv_path)
    
    # Extract base filename for output naming
    base_name = Path(analysis_csv_path).stem.replace('_analysis', '')
    
    # Extract coordinates
    nuclei_centers = []
    for _, row in data_3d.iterrows():
        # Using consistent column names from analyze_masks function
        x_coord = row['x_pos']
        y_coord = row['y_pos'] 
        z_coord = row['z_pos']
        nuclei_centers.append([x_coord, y_coord, z_coord])
    
    nuclei_centers = np.array(nuclei_centers)
    
    # Standardize coordinates
    scaler = StandardScaler()
    nuclei_centers = scaler.fit_transform(nuclei_centers)
    
    # Scale Z axis to account for typical differences in Z resolution
    nuclei_centers[:, 2] *= z_scale_factor
    
    # Perform HDBSCAN clustering
    clusterer = HDBSCAN(min_cluster_size=10, min_samples=8, 
                        cluster_selection_epsilon=6, 
                        cluster_selection_method="eom")
    labels = clusterer.fit_predict(nuclei_centers)
    
    # Count clusters and nuclei
    num_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(f"Predicted number of clusters: {num_clusters}, Total nuclei: {len(nuclei_centers)}")
    
    # Count nuclei in each cluster
    label_counts = Counter(labels)
    total_nuclei = len(nuclei_centers)
    
    # Calculate unclustered nuclei percentage
    unclustered_count = label_counts.get(-1, 0)
    unclustered_percentage = (unclustered_count / total_nuclei) * 100 if total_nuclei > 0 else 0
    
    # Report nuclei per layer
    print("Nuclei per layer:")
    for label, count in sorted(label_counts.items()):
        if label != -1:  # Skip unclustered nuclei in this report
            percentage = (count / total_nuclei) * 100
            print(f"Layer {label+1}: {count} nuclei ({percentage:.2f}%)")
    
    print(f"Unclustered nuclei: {unclustered_count} ({unclustered_percentage:.2f}%)")
    
    # Create summary dictionary for this image
    image_summary = {
        'Image': base_name
    }
    
    # Add counts for each layer
    for label in sorted(label for label in set(labels) if label != -1):
        image_summary[f'Layer {label+1}'] = label_counts[label]
    
    # Add unclustered count
    image_summary['Unclustered'] = label_counts.get(-1, 0)
    
    # Add to summary data if provided
    if summary_data is not None:
        summary_data.append(image_summary)
    
    # Unique labels excluding unclustered
    unique_labels = sorted(label for label in set(labels) if label != -1)
    
    # Assign colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_labels)))
    label_color_map = dict(zip(unique_labels, colors))
    
    # Plot 1: XZ Projection
    plt.figure(figsize=(13, 1.75))
    
    # Plot clusters
    for label in unique_labels:
        cluster_points = nuclei_centers[labels == label]
        plt.scatter(
            cluster_points[:, 0], 
            cluster_points[:, 2],
            c=[label_color_map[label]],
            label=f'Layer {label+1}: {label_counts[label]} nuclei ({(label_counts[label] / total_nuclei) * 100:.1f}%)',
            marker='o'
        )
    
    # Add unclustered points if present
    if -1 in label_counts:
        unclustered_points = nuclei_centers[labels == -1]
        plt.scatter(
            unclustered_points[:, 0], 
            unclustered_points[:, 2], 
            c='gray', 
            marker='x',
            label=f'Unclustered nuclei: {unclustered_count} ({unclustered_percentage:.1f}%)'
        )
    
    plt.xlabel('X Coordinate')
    plt.ylabel('Z Coordinate')
    plt.title(f'HDBSCAN Clustering: XZ Projection, {num_clusters} Layers, {total_nuclei} Nuclei')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    
    # Save XZ plot
    xz_plot_path = os.path.join(output_dir, f"{base_name}_cluster_XZ.png")
    plt.savefig(xz_plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved XZ plot to: {xz_plot_path}")
    
    # Plot 2: XY Projection
    plt.figure(figsize=(12, 6))
    
    # Plot clusters in XY projection
    for label in unique_labels:
        cluster_points = nuclei_centers[labels == label]
        plt.scatter(
            cluster_points[:, 0], 
            cluster_points[:, 1],
            c=[label_color_map[label]],
            label=f'Layer {label+1}: {label_counts[label]} nuclei ({(label_counts[label] / total_nuclei) * 100:.1f}%)',
            marker='o'
        )
    
    # Add unclustered points if present
    if -1 in label_counts:
        unclustered_points = nuclei_centers[labels == -1]
        plt.scatter(
            unclustered_points[:, 0], 
            unclustered_points[:, 1], 
            c='gray', 
            marker='x',
            label=f'Unclustered nuclei: {unclustered_count} ({unclustered_percentage:.1f}%)'
        )
    
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'HDBSCAN Clustering: XY Projection, {num_clusters} Layers, {total_nuclei} Nuclei')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    
    # Save XY plot
    xy_plot_path = os.path.join(output_dir, f"{base_name}_cluster_XY.png")
    plt.savefig(xy_plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved XY plot to: {xy_plot_path}")
    
    # Save clustering results as CSV
    cluster_results = pd.DataFrame({
        'nucleus_id': data_3d['nucleus_id'],
        'x_pos': data_3d['x_pos'],
        'y_pos': data_3d['y_pos'],
        'z_pos': data_3d['z_pos'],
        'volume_um3': data_3d['volume_um3'],
        'cluster_label': labels
    })
    
    cluster_csv_path = os.path.join(output_dir, f"{base_name}_clusters.csv")
    cluster_results.to_csv(cluster_csv_path, index=False)
    print(f"Saved clustering results to: {cluster_csv_path}")
    
    return cluster_results, image_summary

def main():
    try:
        print("\n=== 3D Nuclei Analysis Pipeline ===")
        print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # Parse command line arguments
        args = parse_arguments()
        config_path = args.input
        print(f"Using configuration file: {config_path}")
        
        config = load_config(config_path)
        folders = config['folders']
        nuclei_channel = config['nuclei_channel']
        z_scale_factor = config['z_scale_factor']
        
        print(f"Loaded {len(folders)} folders.")
        print(f"Using nuclei channel: {nuclei_channel}")

        # Initialize ImageJ (headless)
        ij = imagej.init('sc.fiji:fiji', mode='headless')
        IJ = jimport('ij.IJ')

        # Get model path from config
        model_basedir = os.path.dirname(config['model_path'])

        # Load StarDist model
        model = StarDist3D(None, name="stardist-1024-v3", basedir=model_basedir)

        for folder in folders:
            processed_dir = os.path.join(folder, "processed")
            os.makedirs(processed_dir, exist_ok=True)
            # Process ND2 files
            extract_and_resize_nuclei(
                folder,
                processed_dir,
                nuclei_channel=nuclei_channel,
                gaussian_sigma=config['gaussian_sigma'],
             mean_radius=config['mean_radius']
            )

        # Run StarDist segmentation
        run_stardist_segmentation(folders, model)

        # Dictionary to collect summary data across all images
        all_summaries = []

        # Analyze segmentation results and perform clustering
        for folder in folders:
            analysis_dir = os.path.join(folder, "analysis")
            clustering_dir = os.path.join(folder, "clustering")
            os.makedirs(analysis_dir, exist_ok=True)
            os.makedirs(clustering_dir, exist_ok=True)
            masks_dir = os.path.join(folder, "masks")

            for file in os.listdir(masks_dir):
                if file.endswith('_mask.tif') and not file.startswith(('.', '_')):
                    mask_path = os.path.join(masks_dir, file)
                    csv_name = file.replace('_mask.tif', '_analysis.csv')
                    csv_path = os.path.join(analysis_dir, csv_name)
                    print(f"Analyzing {file}")
                    df = analyze_masks(mask_path, csv_path, min_volume=0.5, original_size=(1024, 1024))
                    print(f"{len(df)} nuclei found.")
                    
                    # Run the clustering analysis on each analysis file
                    if len(df) > 10:  # Only cluster if we have enough nuclei
                        try:
                            print(f"Running nuclei clustering on {csv_path}")
                            cluster_results, image_summary = cluster_nuclei(csv_path, clustering_dir, 
                                                                       z_scale_factor=z_scale_factor,
                                                                       summary_data=all_summaries)
                            print(f"Clustering complete for {csv_name}")
                        except Exception as e:
                            print(f"Error during clustering of {csv_name}: {e}")
                            logging.exception(f"Clustering failed for {csv_name}")
                    else:
                        print(f"Skipping clustering for {csv_name} - not enough nuclei detected.")
                        # Add empty summary for skipped images
                        base_name = Path(file).stem.replace('_mask', '')
                        all_summaries.append({
                            'Image': base_name,
                            'Unclustered': len(df)
                        })

            # Create and save summary table for the folder
            if all_summaries:
                # Create DataFrame from summaries
                summary_df = pd.DataFrame(all_summaries)
                
                # Ensure all layer columns exist (fill with 0 if missing)
                max_layer = 0
                for summary in all_summaries:
                    for key in summary.keys():
                        if key.startswith('Layer '):
                            try:
                                layer_num = int(key.split(' ')[1])
                                max_layer = max(max_layer, layer_num)
                            except ValueError:
                                pass
                
                # Make sure all layer columns exist
                for i in range(1, max_layer + 1):
                    layer_col = f'Layer {i}'
                    if layer_col not in summary_df.columns:
                        summary_df[layer_col] = 0
                
                # Ensure 'Unclustered' column exists
                if 'Unclustered' not in summary_df.columns:
                    summary_df['Unclustered'] = 0
                
                # Fill NaN values with 0
                summary_df.fillna(0, inplace=True)
                
                # Convert numeric columns to integers
                for col in summary_df.columns:
                    if col != 'Image':
                        summary_df[col] = summary_df[col].astype(int)
                
                # Reorder columns to put Image first, then layers, then unclustered
                cols = ['Image'] + [f'Layer {i}' for i in range(1, max_layer + 1)] + ['Unclustered']
                summary_df = summary_df[cols]
                
                # Save to CSV
                summary_csv_path = os.path.join(folder, "nuclei_layer_summary.csv")
                summary_df.to_csv(summary_csv_path, index=False)
                print(f"\nNuclei layer summary saved to: {summary_csv_path}")

        print("\n=== Processing Complete ===")
        print(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    except Exception as e:
        print(f"Error: {e}")
        logging.exception("Pipeline execution failed")
    finally:
        ij.dispose() if 'ij' in locals() else None
        print("ImageJ instance closed.")

if __name__ == '__main__':
    main()
