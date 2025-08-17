#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import json
import logging
import os
import shutil
import signal
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import imagej
import jpype
import matplotlib
import numpy as np
import scyjava

matplotlib.rcParams["image.interpolation"] = "none"
import matplotlib.pyplot as plt
import pandas as pd
from csbdeep.utils import normalize
from skimage.measure import regionprops_table
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import StandardScaler
from stardist.models import StarDist3D
from tifffile import imread, imwrite

# ----- logging ---------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("nuclei_analysis.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
np.random.seed(6)


# ----- CLI -------------------------------------------------------------------
def parse_args():
    """Parse command line arguments."""
    ap = argparse.ArgumentParser(description="3D Nuclei Analysis Pipeline")
    ap.add_argument(
        "-i", "--input", default="input_paths.json",
        help="Path to JSON config file"
    )
    # Add hidden argument for Fiji-only mode
    ap.add_argument("--fiji-only", help=argparse.SUPPRESS)
    return ap.parse_args()


# ----- configuration ---------------------------------------------------------
def load_config(path: str) -> dict:
    """Load configuration from JSON file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)

    folders = [p for p in raw.get("folder_paths", []) if os.path.isdir(p)]
    if not folders:
        raise ValueError("'folder_paths' must list existing directories")

    return {
        "folders": folders,
        "nuclei_channel": int(raw.get("nuclei_channel", 1)),
        "gaussian_sigma": float(raw.get("gaussian_sigma", 4.0)),
        "mean_radius": int(raw.get("mean_radius", 3)),
        "z_scale_factor": int(raw.get("z_scale_factor", 32)),
        "n_tiles": raw.get("n_tiles", [1, 1, 1]),
        "model_path": raw.get("model_path"),
    }


# ----- execution-mode selector -----------------------------------------------
def choose_execution_mode(folders):
    """Return user-chosen mode (1,2,3). Defaults to 1 if no outputs exist."""
    outputs_exist = any(
        any((Path(f) / d).exists()
        for d in ("processed", "masks", "analysis", "clustering"))
        for f in folders
    )
    if not outputs_exist:
        return 1

    print("\nDetected existing output folders. Choose execution mode:")
    print("  1 – Full analysis from scratch (clears previous results)")
    print("  2 – Start with StarDist segmentation (skip Fiji pre-processing)")
    print("  3 – Post-processing only (quantification + clustering)\n")
    while True:
        choice = input("Enter 1, 2 or 3 ➔ ").strip()
        if choice in ("1", "2", "3"):
            return int(choice)
        print("Invalid – please enter 1, 2 or 3.")


# ----- Fiji channel extraction -----------------------------------------------
def extract_channel_and_filter(
    in_dir: str, out_dir: str, ch: int, gsig: float, mr: int
):
    """Headless ImageJ: extract nuclei channel & apply 3D filters."""

    ij = imagej.init("sc.fiji:fiji", mode="headless")
    IJ = scyjava.jimport("ij.IJ")
    WM = scyjava.jimport("ij.WindowManager")

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    for fname in sorted(os.listdir(in_dir)):
        if (fname.startswith(('.', '_')) or
                not fname.lower().endswith(('.nd2', '.tif', '.tiff'))):
            continue
        src = os.path.join(in_dir, fname)
        dst = os.path.join(out_dir, f"{Path(fname).stem}_nuclei.tif")
        try:
            imp = IJ.openImage(src)
            if imp is None:
                logging.warning(f"Unreadable: {src}")
                continue
            if imp.getNChannels() < ch:
                logging.error(f"{fname}: channel {ch} out of range")
                imp.close()
                continue

            # Duplicate specified channel
            IJ.run(imp, "Duplicate...", f"title=temp duplicate channels={ch}")
            imp_ch = IJ.getImage()
            imp.close()

            IJ.run(imp_ch, "8-bit", "")
            IJ.run(imp_ch, "Gaussian Blur 3D...", f"x={gsig} y={gsig} z={gsig}")
            IJ.run(imp_ch, "Mean 3D...", f"x={mr} y={mr} z={mr}")
            IJ.saveAs(imp_ch, "Tiff", dst)
            imp_ch.close()
            logging.info(f"Saved {dst}")

        except Exception:
            logging.exception(f"Error on {fname}")
            # Close any orphaned ImageJ windows
            for wid in WM.getIDList() or []:
                win = WM.getImage(wid)
                if win:
                    win.close()

    # Tidy up JVM to avoid lingering threads
    try:
        jpype.shutdownJVM()
    except Exception:
        pass


# ----- StarDist segmentation -------------------------------------------------
def stardist_segment(img, model: StarDist3D, n_tiles):
    """Perform StarDist segmentation on image."""
    labels, _ = model.predict_instances(img, n_tiles=n_tiles)
    return labels.astype(np.uint16)


def run_segmentation(folder: str, model: StarDist3D, n_tiles):
    """Run segmentation on all images in folder."""
    proc = Path(folder, "processed")
    masks = Path(folder, "masks")
    masks.mkdir(exist_ok=True)
    
    # Convert n_tiles to tuple if it's a list
    if isinstance(n_tiles, list):
        n_tiles = tuple(n_tiles)
    
    logging.info(f"Starting segmentation with n_tiles: {n_tiles}")
    
    for tif in sorted(
        p for p in proc.iterdir() if p.suffix.lower() in (".tif", ".tiff")
    ):
        try:
            logging.info(f"Processing {tif.name}")
            img = imread(tif)
            img_norm = normalize(img, 1, 99.8, axis=None)
            
            if img_norm.ndim == 2:
                img_norm = img_norm[np.newaxis, ...]
            
            logging.info(f"Image shape: {img_norm.shape}")
            labels = stardist_segment(img_norm, model, n_tiles)
            base = tif.stem
            mask_path = masks / f"{base}_mask.tif"
            imwrite(mask_path, labels)
            logging.info(f"Saved mask: {mask_path} with {labels.max()} nuclei")

            # Quick QC: mid-Z slice with overlay
            mid = labels.shape[0] // 2
            plt.figure(figsize=(10, 5))
            plt.subplot(121)
            plt.imshow(img_norm[mid], cmap="gray")
            plt.axis("off")
            plt.subplot(122)
            plt.imshow(img_norm[mid], cmap="gray")
            plt.imshow(labels[mid], cmap="jet", alpha=0.45)
            plt.axis("off")
            plt.suptitle(f"{base}: {labels.max()} nuclei")
            plt.tight_layout()
            plt.savefig(masks / f"{base}_QC.png", dpi=120, bbox_inches="tight")
            plt.close()
            
        except Exception as e:
            logging.exception(f"Segmentation failed for {tif.name}: {str(e)}")


# ----- quantification + QC projections ---------------------------------------
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


def process_all_images(
    input_dir: Path, output_dir: Path, z_scale: int = 10, point_size: int = 60
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
    point_size: int = 60
) -> dict:
    """
    Cluster nuclei with 4 projections:
    1. Top-left: XZ (clustered)
    2. Top-right: XY (colors from XZ clustering)
    3. Bottom-left: YZ (clustered)
    4. Bottom-right: XY (colors from YZ clustering)
    Returns image statistics dictionary.
    """
    # Load and validate data
    df = pd.read_csv(csv_path)
    if df.empty:
        logging.warning(f"{csv_path.name}: empty file")
        return
    
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
    xz_clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=epsilon,
        n_jobs=-1
    ).fit(xz_coords)
    df["cluster_xz"] = xz_clusterer.labels_
    
    # YZ clustering (Bottom-left)
    yz_coords = coords_scaled[:, [1, 2]]  # Y and Z
    yz_clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_epsilon=epsilon,
        n_jobs=-1
    ).fit(yz_coords)
    df["cluster_yz"] = yz_clusterer.labels_

    xz_counts = Counter(xz_clusterer.labels_)
    yz_counts = Counter(yz_clusterer.labels_)

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
    
    if summary is not None:
        summary.append(image_stats)

    # Create figure
    fig = plt.figure(figsize=(20, 20))
    
    # Color maps
    xz_clusters = sorted([l for l in set(xz_clusterer.labels_) if l != -1])
    yz_clusters = sorted([l for l in set(yz_clusterer.labels_) if l != -1])
    
    xz_colors = plt.cm.viridis(np.linspace(0, 1, len(xz_clusters)))
    yz_colors = plt.cm.plasma(np.linspace(0, 1, len(yz_clusters)))
    
    # 1. Top-left: XZ projection (clustered)
    ax1 = plt.subplot(2, 2, 1)
    for cluster, color in zip(xz_clusters, xz_colors):
        mask = df["cluster_xz"] == cluster
        ax1.scatter(
            original_coords.loc[mask, "x"],
            original_coords.loc[mask, "z"],
            c=[color],
            s=point_size,
            label=f'XZ Cluster {cluster}'
        )
    
    # 2. Top-right: XY projection (colors from XZ clustering)
    ax2 = plt.subplot(2, 2, 2)
    for cluster, color in zip(xz_clusters, xz_colors):
        mask = df["cluster_xz"] == cluster
        ax2.scatter(
            original_coords.loc[mask, "x"],
            original_coords.loc[mask, "y"],
            c=[color],
            s=point_size,
            label=f'XZ Cluster {cluster}'
        )
    
    # 3. Bottom-left: YZ projection (clustered)
    ax3 = plt.subplot(2, 2, 3)
    for cluster, color in zip(yz_clusters, yz_colors):
        mask = df["cluster_yz"] == cluster
        ax3.scatter(
            original_coords.loc[mask, "y"],
            original_coords.loc[mask, "z"],
            c=[color],
            s=point_size,
            label=f'YZ Cluster {cluster}'
        )
    
    # 4. Bottom-right: XY projection (colors from YZ clustering)
    ax4 = plt.subplot(2, 2, 4)
    for cluster, color in zip(yz_clusters, yz_colors):
        mask = df["cluster_yz"] == cluster
        ax4.scatter(
            original_coords.loc[mask, "x"],
            original_coords.loc[mask, "y"],
            c=[color],
            s=point_size,
            label=f'YZ Cluster {cluster}'
        )

    plt.suptitle(
        f"Nuclei Clustering: {csv_path.stem}\n"
        f"XZ Clusters: {len(xz_clusters)}, YZ Clusters: {len(yz_clusters)}",
        fontsize=16
    )
    plt.tight_layout()

    # Save results
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"{csv_path.stem}_clusters.csv", index=False)
    plt.savefig(
        out_dir / f"{csv_path.stem}_cluster_4proj.png", 
        dpi=150, 
        bbox_inches="tight"
    )
    plt.close()

    return image_stats


# ----- Fiji-only processing mode ---------------------------------------------
def run_fiji_processing(folder: str, cfg: dict):
    """Run Fiji processing in isolation with forced resource cleanup."""
    try:
        logging.info(f"Starting Fiji processing for {folder}")
        extract_channel_and_filter(
            folder, 
            os.path.join(folder, "processed"),
            cfg["nuclei_channel"],
            cfg["gaussian_sigma"],
            cfg["mean_radius"]
        )
    finally:
        # Force cleanup of all resources
        import gc
        gc.collect()
        
        # Special handling for JVM
        try:
            import jpype
            if jpype.isJVMStarted():
                jpype.shutdownJVM()
        except Exception:
            pass
        
        # Cleanup multiprocessing resources
        try:
            import multiprocessing
            multiprocessing.active_children()
        except Exception:
            pass
    print("Fiji subprocess about to exit!", flush=True)
    sys.exit(0)


# ----- emergency shutdown ----------------------------------------------------
def force_terminate():
    # Katya What is that?
    """Forcefully terminate processes and exit."""
    try:
        import jpype
        if jpype.isJVMStarted():
            jpype.shutdownJVM()
    except Exception:
        pass
    try:
        import tensorflow as tf
        tf.keras.backend.clear_session()
    except Exception:
        pass
    os._exit(0)


def main():
    """
    3D Nuclei Analysis Pipeline – native resolution (full)   Updated: 29 Jun 2025
    ------------------------------------------------------------------------------
     • Three execution modes (auto-detected if previous results exist)
          1) Full analysis from scratch
          2) StarDist segmentation only (skip Fiji pre-processing)
          3) Post-processing only (quantification + clustering)
     • Fiji-based nuclei-channel extraction & 3D filtering (no resizing)
     • Tile-aware StarDist 3D segmentation (multi-CPU)
     • Per-image quantification → CSV
     • QC plots (mid-Z mask overlay + 3-view colour projections & centroids)
     • HDBSCAN 3-D clustering into nuclear "layers" → per-image + summary CSV
     • Robust cleanup of ImageJ JVM & TensorFlow sessions (script always exits)
    """
    logging.info("=== 3D Nuclei Analysis Pipeline ===")
    logging.info(f"Start: {datetime.now():%Y-%m-%d %H:%M:%S}")
    args = parse_args()

    # Загружаем конфиг
    cfg = load_config(args.input)
    mode = choose_execution_mode(cfg["folders"])
    logging.info(f"Execution mode: { {1:'full',2:'stardist',3:'post'}[mode] }")

    signal.signal(signal.SIGTERM, lambda *a: force_terminate())

    try:
        summaries = []
        model = None

        # Загружаем модель для StarDist, если нужна
        if mode in (1, 2):
            if cfg["model_path"]:
                mdir = Path(cfg["model_path"]).resolve()
                model = StarDist3D(None, name=mdir.name, basedir=str(mdir.parent))
                logging.info(f"Loaded custom model: {mdir}")
            else:
                model = StarDist3D.from_pretrained("3D_demo")
                logging.info("Using built-in 3D_demo model")

        for fld in cfg["folders"]:
            fld = Path(fld)
            proc = fld / "processed"
            masks = fld / "masks"
            anal = fld / "analysis"
            clus = fld / "clustering"

            # Режим 1: с нуля — чистим папки и запускаем Fiji/ImageJ
            if mode == 1:
                for d in (proc, masks, anal, clus):
                    if d.exists():
                        shutil.rmtree(d)
                for d in (proc, masks, anal, clus):
                    d.mkdir(parents=True, exist_ok=True)
                logging.info(f"Running Fiji processing in main process for {fld}")
                run_fiji_processing(str(fld), cfg)
                logging.info(f"Fiji processing completed for {fld}")

            # В режимах 1 и 2 — запускаем StarDist
            if mode in (1, 2):
                logging.info(f"Starting segmentation with n_tiles: {cfg['n_tiles']}")
                run_segmentation(str(fld), model, cfg["n_tiles"])

            # Квантификация и кластеризация
            for msk in masks.glob("*_mask.tif"):
                csv_out = anal / msk.name.replace("_mask.tif", "_analysis.csv")
                qc_png = anal / msk.name.replace("_mask.tif", "_3D_QC.png")
                if not csv_out.exists() or mode in (1, 2):
                    quantify(msk, csv_out, qc_png)
                cluster(csv_out, clus, cfg["z_scale_factor"], summaries)

            # Финальная обработка всех изображений в папке
            process_all_images(clus, clus, cfg["z_scale_factor"])

        # Итоговая сводка по всем папкам
        if summaries:
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

    except Exception as e:
        logging.exception(f"Pipeline failed: {e}")
    finally:
        force_terminate()


# Handle Fiji-only subprocess mode
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logging.warning("Interrupted by user")
        sys.exit(1)
    except Exception:
        logging.exception("Fatal error")
        sys.exit(1)