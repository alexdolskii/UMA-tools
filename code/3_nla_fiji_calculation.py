#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fiji Calculation Module
Performs quantification, QC projections, and clustering analysis
"""

import argparse
import logging
from collections import Counter
from pathlib import Path
from datetime import datetime
from nla_config_loader import load_config
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import StandardScaler
from tifffile import imread


def generate_qc_projections(mask_arr, df, out_png):
    """
    Create color projections (XY,XZ,YZ) and centroid scatter plots for QC.
    """
    lbl_max = mask_arr.max()
    if lbl_max == 0:
        return
    rng = np.random.RandomState(42)
    cmap = rng.rand(lbl_max + 1, 3)
    cmap[0] = 0  # Background → black

    z, y, x = mask_arr.shape
    rgb = np.zeros((z, y, x, 3), dtype=np.float32)
    for lbl in range(1, lbl_max + 1):
        rgb[mask_arr == lbl] = cmap[lbl]

    xy = rgb.max(axis=0)
    xz = np.flipud(rgb.max(axis=1))
    yz = np.flipud(rgb.max(axis=2))

    plt.figure(figsize=(18, 10))
    # XY projection
    plt.subplot(231)
    plt.imshow(xy)
    plt.title("XY projection")
    plt.axis("off")
    # XZ projection
    plt.subplot(232)
    plt.imshow(xz)
    plt.title("XZ projection")
    plt.axis("off")
    # YZ projection
    plt.subplot(233)
    plt.imshow(yz)
    plt.title("YZ projection")
    plt.axis("off")
    # Centroids XY
    plt.subplot(234)
    plt.imshow(xy)
    plt.axis("off")
    plt.scatter(
        df["x"], df["y"], s=12, c="white", linewidths=0.3, edgecolors="k"
    )
    plt.title("Centroids (XY)")
    # Centroids XZ
    plt.subplot(235)
    plt.imshow(xz)
    plt.axis("off")
    plt.scatter(
        df["x"], z - df["z"], s=12, c="white", linewidths=0.3, edgecolors="k"
    )
    plt.title("Centroids (XZ)")
    # Centroids YZ
    plt.subplot(236)
    plt.imshow(yz)
    plt.axis("off")
    plt.scatter(
        df["y"], z - df["z"], s=12, c="white", linewidths=0.3, edgecolors="k"
    )
    plt.title("Centroids (YZ)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close()


def quantify(
    mask_path: Path, csv_out: Path, qc_png: Path, min_vox: int = 10
) -> pd.DataFrame:
    """Quantify nuclei properties and generate CSV + QC plot."""
    labels = imread(mask_path)
    props = regionprops_table(
        labels,
        properties=("label", "area", "centroid", "equivalent_diameter_area")
    )
    df = (pd.DataFrame(props)
          .rename(columns={
              "label": "nucleus_id",
              "area": "vol_vox",
              "centroid-0": "z",
              "centroid-1": "y",
              "centroid-2": "x",
              "equivalent_diameter_area": "diam_px"
          }))
    df = df[df["vol_vox"] >= min_vox]
    df.to_csv(csv_out, index=False)

    # Generate QC plots
    try:
        generate_qc_projections(labels, df, qc_png)
    except Exception:
        logging.warning(f"QC plot failed for {mask_path.name}")

    return df


def cluster(
        csv_path: Path,
        out_dir: Path,
        z_scale: int = 10,
        summary: list = None,
        *,
        trim_method: str = "MAD",
        k: float | tuple = 2.5,
        min_cluster_size: int = 5,
        min_samples: int = 3,
        epsilon: float = 4.0,
        min_nuclei: int = 10,
        normalize: bool = True,
        point_size: int = 60) -> dict:
    """
    Cluster nuclei with 4 projections:
    1. Top-left: XZ (clustered)
    2. Top-right: XY (colors from XZ clustering)
    3. Bottom-left: YZ (clustered)
    4. Bottom-right: XY (colors from YZ clustering)
    Returns image statistics dictionary.
    Always plot all points; unclustered (-1) are gray.
    """
    # Load and validate data
    df = pd.read_csv(csv_path)
    if df.empty:
        logging.warning(f"{csv_path.name}: empty file")
        return
        
    if len(df) < 3:
        logging.warning(f"{csv_path.name}: Only {len(df)} nuclei found. Clustering requires at least 3 nuclei.")
        # Still proceed but will skip clustering

    if summary is None:
        summary = []

    # Store original coordinates for visualization
    original_coords = df[["x", "y", "z"]].copy()

    # Prepare scaled coordinates for clustering
    coords = df[["x", "y", "z"]].to_numpy()
    if normalize:
        coords = StandardScaler().fit_transform(coords)
    coords_scaled = coords.copy()
    coords_scaled[:, 2] *= z_scale  # Scale Z-axis only

    # XZ clustering (Top-left)
    xz_coords = coords_scaled[:, [0, 2]]  # X and Z
    
    # Adjust clustering parameters based on data size
    actual_min_samples = min(min_samples, len(xz_coords))
    actual_min_cluster_size = min(min_cluster_size, len(xz_coords))
    
    if len(xz_coords) < 3:
        # Not enough data for clustering, assign all to noise
        logging.warning(f"{csv_path.name}: Only {len(xz_coords)} nuclei found, skipping XZ clustering")
        df["cluster_xz"] = [-1] * len(df)
        xz_clusterer = None
    else:
        xz_clusterer = HDBSCAN(
            min_cluster_size=actual_min_cluster_size,
            min_samples=actual_min_samples,
            cluster_selection_epsilon=epsilon,
            n_jobs=-1
        ).fit(xz_coords)
        df["cluster_xz"] = xz_clusterer.labels_
    
    # YZ clustering (Bottom-left)
    yz_coords = coords_scaled[:, [1, 2]]  # Y and Z
    
    if len(yz_coords) < 3:
        # Not enough data for clustering, assign all to noise
        logging.warning(f"{csv_path.name}: Only {len(yz_coords)} nuclei found, skipping YZ clustering")
        df["cluster_yz"] = [-1] * len(df)
        yz_clusterer = None
    else:
        yz_clusterer = HDBSCAN(
            min_cluster_size=actual_min_cluster_size,
            min_samples=actual_min_samples,
            cluster_selection_epsilon=epsilon,
            n_jobs=-1
        ).fit(yz_coords)
        df["cluster_yz"] = yz_clusterer.labels_

    # Handle cases where clustering wasn't performed
    if xz_clusterer is not None:
        xz_counts = Counter(xz_clusterer.labels_)
    else:
        xz_counts = Counter([-1] * len(df))
        
    if yz_clusterer is not None:
        yz_counts = Counter(yz_clusterer.labels_)
    else:
        yz_counts = Counter([-1] * len(df))

    # Prepare statistics for summary
    image_stats = {
        'Image': csv_path.stem,
        'Total_Nuclei': len(df),
        'XZ_Clusters': len([k for k in xz_counts if k != -1]),
        'YZ_Clusters': len([k for k in yz_counts if k != -1]),
        'XZ_Unclustered': xz_counts.get(-1, 0),
        'YZ_Unclustered': yz_counts.get(-1, 0),
        'XZ_Cluster_Counts': dict(xz_counts),
        'YZ_Cluster_Counts': dict(yz_counts)
    }
    summary.append(image_stats)

    # Create figure
    fig = plt.figure(figsize=(20, 20))

    # Color maps
    if xz_clusterer is not None:
        xz_clusters = sorted([l for l in set(xz_clusterer.labels_) if l != -1])
    else:
        xz_clusters = []
        
    if yz_clusterer is not None:
        yz_clusters = sorted([l for l in set(yz_clusterer.labels_) if l != -1])
    else:
        yz_clusters = []

    xz_colors = plt.cm.viridis(np.linspace(0, 1, max(1, len(xz_clusters))))[:len(xz_clusters)]
    yz_colors = plt.cm.plasma(np.linspace(0, 1, max(1, len(yz_clusters))))[:len(yz_clusters)]

    # masks for noise
    xz_noise = (df["cluster_xz"] == -1)
    yz_noise = (df["cluster_yz"] == -1)

    # 1. Top-left: XZ projection (clustered)
    ax1 = plt.subplot(2, 2, 1)
    if xz_noise.any():
        ax1.scatter(
            original_coords.loc[xz_noise, "x"],
            original_coords.loc[xz_noise, "z"],
            c="0.8", s=point_size, label="Unclustered", zorder=1
        )
    for cluster, color in zip(xz_clusters, xz_colors):
        m = (df["cluster_xz"] == cluster)
        ax1.scatter(
            original_coords.loc[m, "x"],
            original_coords.loc[m, "z"],
            c=[color], s=point_size, label=f'XZ Cluster {cluster}', zorder=2
        )
    ax1.set_title("XZ (HDBSCAN)")
    ax1.set_xlabel("X"); ax1.set_ylabel("Z")

    # 2. Top-right: XY projection (colors from XZ clustering)
    ax2 = plt.subplot(2, 2, 2)
    if xz_noise.any():
        ax2.scatter(
            original_coords.loc[xz_noise, "x"],
            original_coords.loc[xz_noise, "y"],
            c="0.8", s=point_size, label="Unclustered", zorder=1
        )
    for cluster, color in zip(xz_clusters, xz_colors):
        m = (df["cluster_xz"] == cluster)
        ax2.scatter(
            original_coords.loc[m, "x"],
            original_coords.loc[m, "y"],
            c=[color], s=point_size, label=f'XZ Cluster {cluster}', zorder=2
        )
    ax2.set_title("XY (colored by XZ)")
    ax2.set_xlabel("X"); ax2.set_ylabel("Y")

    # 3. Bottom-left: YZ projection (clustered)
    ax3 = plt.subplot(2, 2, 3)
    if yz_noise.any():
        ax3.scatter(
            original_coords.loc[yz_noise, "y"],
            original_coords.loc[yz_noise, "z"],
            c="0.8", s=point_size, label="Unclustered", zorder=1
        )
    for cluster, color in zip(yz_clusters, yz_colors):
        m = (df["cluster_yz"] == cluster)
        ax3.scatter(
            original_coords.loc[m, "y"],
            original_coords.loc[m, "z"],
            c=[color], s=point_size, label=f'YZ Cluster {cluster}', zorder=2
        )
    ax3.set_title("YZ (HDBSCAN)")
    ax3.set_xlabel("Y"); ax3.set_ylabel("Z")

    # 4. Bottom-right: XY projection (colors from YZ clustering)
    ax4 = plt.subplot(2, 2, 4)
    if yz_noise.any():
        ax4.scatter(
            original_coords.loc[yz_noise, "x"],
            original_coords.loc[yz_noise, "y"],
            c="0.8", s=point_size, label="Unclustered", zorder=1
        )
    for cluster, color in zip(yz_clusters, yz_colors):
        m = (df["cluster_yz"] == cluster)
        ax4.scatter(
            original_coords.loc[m, "x"],
            original_coords.loc[m, "y"],
            c=[color], s=point_size, label=f'YZ Cluster {cluster}', zorder=2
        )
    ax4.set_title("XY (colored by YZ)")
    ax4.set_xlabel("X"); ax4.set_ylabel("Y")

    plt.suptitle(
        f"Nuclei Clustering: {csv_path.stem}\n"
        f"XZ Clusters: {len(xz_clusters)}, YZ Clusters: {len(yz_clusters)}",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{csv_path.stem}_clusters.csv", index=False)
    plt.savefig(out_dir / f"{csv_path.stem}_cluster_4proj.png", dpi=150, bbox_inches="tight")
    plt.close()

    return image_stats


def create_summary_dataframe(summary_data: list) -> pd.DataFrame:
    """Create detailed summary DataFrame from image processing results."""
    rows = []

    for image_data in summary_data:
        row = {
            'Image': image_data['Image'],
            'Total_Nuclei': image_data['Total_Nuclei'],
            'XZ_Clusters': image_data.get('XZ_Clusters', 0),
            'YZ_Clusters': image_data.get('YZ_Clusters', 0),
            'XZ_Unclustered': image_data.get('XZ_Unclustered', 0),
            'YZ_Unclustered': image_data.get('YZ_Unclustered', 0)
        }

        # Add cluster counts for XZ
        xz_counts = image_data.get('XZ_Cluster_Counts', {})
        for cluster, count in xz_counts.items():
            if cluster != -1:
                row[f'XZ_Cluster_{cluster}_Count'] = count

        # Add cluster counts for YZ
        yz_counts = image_data.get('YZ_Cluster_Counts', {})
        for cluster, count in yz_counts.items():
            if cluster != -1:
                row[f'YZ_Cluster_{cluster}_Count'] = count

        rows.append(row)

    return pd.DataFrame(rows)


def process_all_images(
        input_dir: Path,
        output_dir: Path,
        z_scale: int = 10,
        point_size: int = 60
):
    """
    Process all CSV files in input_dir and create a summary CSV.
    Uses cluster() for processing and create_summary_dataframe() for reporting.
    """
    summary_data = []

    for csv_file in input_dir.glob('*_analysis.csv'):
        # Process each image
        image_stats = cluster(
            csv_path=csv_file,
            out_dir=output_dir,
            z_scale=z_scale,
            summary=summary_data,
            point_size=point_size
        )

    # Create comprehensive summary
    if summary_data:
        summary_df = create_summary_dataframe(summary_data)
        summary_csv = output_dir / 'clustering_summary.csv'
        summary_df.to_csv(summary_csv, index=False)
        logging.info(f"Saved comprehensive summary to {summary_csv}")


def run_quantification_and_clustering(folder_path: str, config: dict):
    """Run quantification and clustering for a single folder."""
    fld = Path(folder_path)
    masks = fld / "masks"
    anal = fld / "analysis"
    clus = fld / "clustering"

    summaries = []

    # Quantification
    for msk in masks.glob("*_mask.tif"):
        csv_out = anal / msk.name.replace("_mask.tif", "_analysis.csv")
        qc_png = anal / msk.name.replace("_mask.tif", "_3D_QC.png")
        quantify(msk, csv_out, qc_png)

    # Clustering
    for csv_file in anal.glob("*_analysis.csv"):
        cluster(csv_file, clus, config["z_scale_factor"], summaries)

    # Process all images and create summary
    process_all_images(clus, clus, config["z_scale_factor"])

    return summaries


def main(input_file_path):
    """
    3D Nuclei Analysis Pipeline
    ------------------------------------------------------------------------------
     The main function to implement the 3rd step of
     the 3D Nuclei Analysis Pipeline.

    The function implements post-processing:
    quantification + clustering
    """
    np.random.seed(6)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s")
    print("=== 3D Nuclei Analysis Pipeline ===")
    print("3 - Post-processing")
    print(f"Start: {datetime.now():%Y-%m-%d %H:%M:%S}")

    cfg = load_config(input_file_path)

    summaries = []

    for fld in cfg["folders"]:
        fld = Path(fld)
        proc = fld / "processed"
        masks = fld / "masks"
        anal = fld / "analysis"
        clus = fld / "clustering"

        fh = logging.FileHandler(os.path.join(fld,
                                              "nuclei_analysis3.log"), mode='w')
        fh.setLevel(logging.INFO)
        logging.getLogger('').addHandler(fh)

        for msk in masks.glob("*_mask.tif"):
            csv_out = anal / msk.name.replace("_mask.tif", "_analysis.csv")
            qc_png = anal / msk.name.replace("_mask.tif", "_3D_QC.png")
            if not csv_out.exists():
                quantify(msk, csv_out, qc_png)
            cluster(csv_out, clus, cfg["z_scale_factor"], summaries)

        # Финальная обработка всех изображений в папке
        process_all_images(clus, clus, cfg["z_scale_factor"])

    # Итоговая сводка по всем папкам
    if len(summaries) > 0:
        logging.info(f"Collected data for {len(summaries)} images")
        summary_df = create_summary_dataframe(summaries)
        if cfg["folders"]:
            first_folder = Path(cfg["folders"][0])
            clus_dir = first_folder / "clustering"
            summary_path = clus_dir / "clustering_summary.csv"
            summary_df.to_csv(summary_path, index=False)
            logging.info(f"Saved final summary to {summary_path}")
        else:
            summary_df.to_csv("clustering_summary.csv", index=False)
            logging.info("Saved final summary to current directory")

    logging.info("=== Processing Complete ===")
    logging.info(f"End: {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3D Nuclei Analysis Pipeline")
    parser.add_argument("-i",
                        "--input",
                        default="input_paths.json",
                        help="Path to JSON config file"
    )
    args = parser.parse_args()
    main(args.input)
