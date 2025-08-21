#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fiji Channel Extraction Module
Extracts nuclei channel and applies 3D filters using headless ImageJ/Fiji
"""

import logging
import os
from pathlib import Path
import argparse
import shutil
from datetime import datetime
import numpy as np

import imagej
import scyjava
from nla_config_loader import load_config


def extract_channel_and_filter(
    in_dir: str,
    out_dir: str,
    ch: int,
    gsig: float,
    mr: int) -> None:
    """
    Headless ImageJ: extract nuclei channel
    and apply 3D filters.
    """
    ij = imagej.init("sc.fiji:fiji", mode="headless")
    IJ = scyjava.jimport("ij.IJ")
    WM = scyjava.jimport("ij.WindowManager")

    try:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        for fname in sorted(os.listdir(in_dir)):
            if (fname.startswith(('.', '_')) or
                    not fname.lower().endswith(('.nd2',
                                                '.tif',
                                                '.tiff'))):
                continue
            src = os.path.join(in_dir, fname)
            dst = os.path.join(out_dir,
                               f"{Path(fname).stem}_nuclei.tif")

            imp = IJ.openImage(src)
            if imp is None:
                logging.warning(f"Unreadable: {src}")
                continue
            if imp.getNChannels() < ch:
                logging.error(f"{fname}: channel {ch} out of range")
                imp.close()
                continue

            # Duplicate specified channel
            IJ.run(imp, "Duplicate...",
                   f"title=temp duplicate channels={ch}")
            imp_ch = IJ.getImage()
            imp.close()

            IJ.run(imp_ch, "8-bit", "")
            IJ.run(imp_ch,
                   "Gaussian Blur 3D...",
                   f"x={gsig} y={gsig} z={gsig}")
            IJ.run(imp_ch, "Mean 3D...",
                   f"x={mr} y={mr} z={mr}")
            IJ.saveAs(imp_ch, "Tiff", dst)
            imp_ch.close()
            logging.info(f"Saved {dst}")
    finally:
        for wid in WM.getIDList() or []:
            win = WM.getImage(wid)
            if win:
                win.close()
        ij.dispose()
        scyjava.jimport('java.lang.System').exit(0)


def main_fiji_processing(input_file_path):
    """
    3D Nuclei Analysis Pipeline
    ------------------------------------------------------------------------------
     The main function to implement the 1st step of
     the 3D Nuclei Analysis Pipeline.

    The function extracts nuclei channel and
    applies 3D filters using headless ImageJ/Fiji
    """
    np.random.seed(6)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s")

    print("=== 3D Nuclei Analysis Pipeline ===")
    print("1 - Fiji pre-processing")
    print(f"Start: {datetime.now():%Y-%m-%d %H:%M:%S}")

    cfg = load_config(input_file_path)

    for fld in cfg["folders"]:
        fld = Path(fld)
        proc = fld / "processed"
        masks = fld / "masks"
        anal = fld / "analysis"
        clus = fld / "clustering"

        fh = logging.FileHandler(os.path.join(fld,
                                              "nuclei_analysis1.log"), mode='w')
        fh.setLevel(logging.INFO)
        logging.getLogger('').addHandler(fh)

        outputs_exist = any(
            any((Path(f) / d).exists()
                for d in ("processed", "masks", "analysis", "clustering"))
            for f in cfg["folders"]
        )
        if outputs_exist:
            start_analysis = (input("\nThe files from the first step exists. "
                                    "Do you want to rewrite them? (y/n): ")
                                    .strip().lower())
            if start_analysis in ('no', 'n'):
                raise ValueError("Analysis canceled by user.")
            elif start_analysis not in ('yes', 'y', 'no', 'n'):
                raise ValueError("Incorrect input. Please enter y/n or yes/no")

        for d in (proc, masks, anal, clus):
            if d.exists():
                shutil.rmtree(d)
        for d in (proc, masks, anal, clus):
            d.mkdir(parents=True, exist_ok=True)
        logging.info(f"Running Fiji processing in main process for {fld}")
        logging.info(f"Starting Fiji processing for {str(fld)}")

        extract_channel_and_filter(
            str(fld),
            os.path.join(str(fld), "processed"),
            cfg["nuclei_channel"],
            cfg["gaussian_sigma"],
            cfg["mean_radius"]
        )
        logging.info(f"Fiji processing completed for {fld}")

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
    main_fiji_processing(args.input)
