#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
StarDist Prediction Module
Performs 3D nuclei segmentation using StarDist models
"""

import argparse
import logging
from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from csbdeep.utils import normalize
from stardist.models import StarDist3D
from tifffile import imread, imwrite
from nla_config_loader import load_config


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


def load_stardist_model(model_path: str = None):
    """Load StarDist model from path or use built-in model."""
    if model_path:
        mdir = Path(model_path).resolve()
        model = StarDist3D(None, name=mdir.name, basedir=str(mdir.parent))
        logging.info(f"Loaded custom model: {mdir}")
    else:
        model = StarDist3D.from_pretrained("3D_demo")
        logging.info("Using built-in 3D_demo model")

    return model


def main_star_dist_segmentation(input_file_path):
    """
    3D Nuclei Analysis Pipeline
    ------------------------------------------------------------------------------
     The main function to implement the 2nd step of
     the 3D Nuclei Analysis Pipeline.

    The function performs tile-aware StarDist
    3D segmentation (multi-CPU or GPU)
    """
    np.random.seed(6)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s")
    print("=== 3D Nuclei Analysis Pipeline ===")
    print("2 - StarDist segmentation")
    print(f"Start: {datetime.now():%Y-%m-%d %H:%M:%S}")

    cfg = load_config(input_file_path)

    model = None
    model = load_stardist_model(cfg["model_path"])

    for fld in cfg["folders"]:
        fld = Path(fld)
        proc = fld / "processed"
        masks = fld / "masks"
        anal = fld / "analysis"
        clus = fld / "clustering"

        fh = logging.FileHandler(os.path.join(fld,
                                              "nuclei_analysis2.log"), mode='w')
        fh.setLevel(logging.INFO)
        logging.getLogger('').addHandler(fh)

        logging.info(f"Starting segmentation with n_tiles: {cfg['n_tiles']}")
        run_segmentation(str(fld), model, cfg["n_tiles"])

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
    main_star_dist_segmentation(args.input)
