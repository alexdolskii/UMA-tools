#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fibronectin‐fiber alignment assay, ImageJ + OrientationJ version
---------------------------------------------------------------
Workflow (unchanged in essence):

1.  Read a JSON file that lists folders to process.
2.  **Part 1** – open every *.tif / .tiff / .nd2 / .oif / .oib* stack,
   extract the fibronectin channel, generate a max-intensity Z-projection,
   resize to a standard size, and save as “…_processed.tif”.
3.  **Part 2** – feed every processed TIFF to the *OrientationJ Analysis*
   and *OrientationJ Distribution* plugins (inside ImageJ) and save:
     • an HSV survey image “…_oj_analysis.tif”
     • the distribution table “…_oj_distribution.csv”
4.  **Part 3** – merge all CSV tables, classify each field of view
   (‘aligned’ vs ‘disorganized’), and write both individual and summary CSVs.

Folder layout created in each input folder:

    Alignment_assay_results_angle_<ANGLE>_<TIMESTAMP>/
        ├─ Images/              (survey images)
        ├─ Tables/              (OrientationJ CSVs)
        ├─ Analysis/            (post-processed CSVs + summary)
        └─ log.log              (per-folder log)

No image-processing logic (OrientationJ calls, etc.) has been modified.
Only *additional* handling for *.oif/.oib* has been inserted.

--------------------------------------------------------------------
"""

import argparse
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple
import sys
import imagej
import pandas as pd
import scyjava as sj


# ──────────────────────────────── Exceptions ────────────────────────────────

class ImageJInitializationError(RuntimeError):
    """Raised when ImageJ fails to initialise."""


# ────────────────────────── Generic helper functions ─────────────────────────

def correct_angle(angle: float) -> float:
    """Map any angle into the range [-90°, +90°]."""
    if angle < -90:
        return angle + 180
    if angle > 90:
        return angle - 180
    return angle


def read_folder_list(json_path: str) -> List[str]:
    """
    Parse *input_paths.json* and sanity-check that the folders exist.
    The file must look like:
        {"folder_paths": ["C:/data/groupA", "C:/data/groupB", ...]}
    """
    if not os.path.isfile(json_path):
        raise FileNotFoundError(f"File '{json_path}' does not exist.")

    with open(json_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    raw_paths = data.get("folder_paths", [])
    if not raw_paths:
        raise ValueError("JSON does not contain key 'folder_paths'.")

    valid_paths: List[str] = []
    for p in raw_paths:
        if os.path.isdir(p):
            n_files = len(os.listdir(p))
            exts = {Path(f).suffix.lower() for f in os.listdir(p)}
            print(f"Folder: {p}\n  files: {n_files}\n  types: {', '.join(exts)}")
            valid_paths.append(p)
        else:
            logging.warning(f"Folder '{p}' does not exist – skipped.")
    if not valid_paths:
        raise ValueError("No existing folders to process.")
    print(f"\nFound {len(valid_paths)} folder(s) for processing.")
    return valid_paths


def make_result_dirs(root: str, angle_str: str) -> Tuple[str, str, str]:
    """Create results/, Tables/, Images/ under *root*; return their paths."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    results = Path(root) / f"Alignment_assay_results_angle_{angle_str}_{ts}"
    tables = results / "Tables"
    images = results / "Images"
    for d in (results, tables, images):
        d.mkdir(parents=True, exist_ok=True)
    # one log per input folder
    logging.FileHandler(results / "log.log", mode="w")
    print(f"Results will be written to: {results}")
    return str(results), str(tables), str(images)


# ───────────────────────────── ImageJ utilities ─────────────────────────────

def init_imagej() -> "imagej.ImageJ":
    # === 1. Determine the location of the script and look for the Fiji folder in the same directory. ===
    SCRIPT_DIR = Path(__file__).resolve().parent          # code/
    FIJI_DIR = (SCRIPT_DIR.parent / "Fiji.app").resolve()

    if not FIJI_DIR.exists():
        raise FileNotFoundError(
            f"Fiji.app not found at {FIJI_DIR}. "
            "Place the portable Fiji.app one level above the 'code' folder "
            "or pass --fiji <path>."
        )   
    # Initialize ImageJ in interactive mode
    print("Initializing ImageJ...")

    ij = imagej.init(str(FIJI_DIR), mode='interactive')

    print("ImageJ initialization completed.")
    return ij


def close_all_images() -> None:
    """Utility to close every open ImageJ window."""
    IJ = sj.jimport("ij.IJ")
    IJ.run("Close All")


# ───────────────────────────────── Part 1 ────────────────────────────────────

def part1_z_projection(
    ij: "imagej.ImageJ",
    folder: str,
    results_dir: str,
    fn_channel: int,
    target_w: int,
    target_h: int,
) -> Dict[str, Dict]:
    """
    Generate resized max-intensity projections of the fibronectin channel.
    Returns a dict with basic Z-stack info keyed by processed filename.
    """
    IJ = sj.jimport("ij.IJ")
    ZProjector = sj.jimport("ij.plugin.ZProjector")

    z_info: Dict[str, Dict] = {}

    for f in os.listdir(folder):
        if not f.lower().endswith((".tif", ".tiff", ".nd2", ".oif", ".oib")):
            continue
        full_path = os.path.join(folder, f)
        print(f"\n[Part 1] Opening: {full_path}")
        close_all_images()
        imp = IJ.openImage(full_path)
        if imp is None:
            logging.warning(f"Could not open '{f}' – skipped.")
            continue

        w, h, c, z, t = imp.getDimensions()
        if fn_channel > c:
            logging.warning(f"{f}: requested channel {fn_channel} > {c}")
            imp.close(); continue

        # Extract fibronectin channel (ImageJ channels are 1-based)
        imp.setC(fn_channel)
        IJ.run(imp, "Duplicate...", f"title=fibro channels={fn_channel}")
        ch_img = IJ.getImage()

        # Z-projection
        zp = ZProjector(ch_img)
        zp.setMethod(ZProjector.MAX_METHOD)
        zp.doProjection()
        proj = zp.getProjection()
        proj = proj.resize(target_w, target_h, "bilinear")
        IJ.run(proj, "8-bit", "")

        out_name = Path(f).stem + "_processed.tif"
        out_path = os.path.join(results_dir, out_name)
        IJ.saveAs(proj, "Tiff", out_path)
        logging.info(f"Saved projection → {out_path}")
        proj.close(); ch_img.close(); imp.close(); close_all_images()

        z_info[Path(out_name).stem] = {
            "original_filename": f,
            "number_of_z_stacks": z,
            "z_stack_type": "slices" if z > 1 else "channels",
        }
    return z_info


# ───────────────────────────────── Part 2 ────────────────────────────────────

def part2_orientationj(results_dir: str, images_dir: str) -> None:
    """
    Run OrientationJ Analysis + Distribution on every *_processed.tif file.
    Saves the HSV survey image and distribution CSV (as in the legacy script).
    """
    IJ = sj.jimport("ij.IJ")
    WM = sj.jimport("ij.WindowManager")

    targets = [
        f for f in os.listdir(results_dir) if f.lower().endswith("_processed.tif")
    ]
    if not targets:
        logging.warning("Part 2: no processed TIFFs found.")
        return

    for fname in targets:
        path = os.path.join(results_dir, fname)
        print(f"\n[Part 2] OrientationJ on: {fname}")
        close_all_images()
        imp = IJ.openImage(path)
        if imp is None:
            logging.warning(f"Could not reopen '{fname}'.")
            continue
        imp.show()

        # Analysis
        IJ.run(
            "OrientationJ Analysis",
            "tensor=3.0 gradient=4 color-survey=on "
            "hsb=on hue=Orientation sat=Coherency bri=Original-Image radian=on",
        )
        IJ.wait(500)
        survey = WM.getImage("OJ-Color-survey-1")
        if survey:
            survey_out = os.path.join(
                images_dir, Path(fname).stem + "_oj_analysis.tif"
            )
            IJ.saveAs(survey, "Tiff", survey_out)
            survey.close()
            logging.info(f"Saved survey image → {survey_out}")

        close_all_images()
        # Distribution
        imp = IJ.openImage(path); imp.show()
        IJ.run(
            "OrientationJ Distribution",
            "tensor=3.0 gradient=4 radian=on histogram=on table=on "
            "min-coherency=0.0 min-energy=0.0",
        )
        IJ.wait(500)
        csv_out = os.path.join(
            results_dir, "Tables", Path(fname).stem + "_oj_distribution.csv"
        )
        try:
            IJ.saveAs("Results", csv_out)
            logging.info(f"Saved distribution table → {csv_out}")
            IJ.run("Clear Results")
        except Exception as e:
            logging.error(f"Could not save results for '{fname}': {e}")
        close_all_images()


# ───────────────────────────────── Part 3 ────────────────────────────────────

def part3_summarise(
    results_dir: str,
    analysis_dir: str,
    angle: float,
    z_info: Dict[str, Dict],
) -> None:
    """Post-process OrientationJ CSVs and produce a summary table."""
    tables_dir = Path(results_dir) / "Tables"
    csvs = list(tables_dir.glob("*.csv"))
    if not csvs:
        logging.warning("Part 3: no CSVs to process.")
        return

    analysis_dir = Path(analysis_dir)
    analysis_dir.mkdir(exist_ok=True)

    summary_rows: List[Dict] = []

    for csv_path in csvs:
        df = pd.read_csv(csv_path)
        df.columns = ["ori_angle", "occ_value"]
        peak_angle = df.loc[df["occ_value"].idxmax(), "ori_angle"]

        df["angles_norm"] = df["ori_angle"] - peak_angle
        df["angles_corr"] = df["angles_norm"].apply(correct_angle)
        df["rank"] = df["angles_corr"].rank(method="min")
        df["perc"] = df["occ_value"] / df["occ_value"].sum() * 100

        mask = df["angles_corr"].between(-angle, angle)
        pct_aligned = df.loc[mask, "perc"].sum()
        mode = "aligned" if pct_aligned >= 55 else "disorganized"

        # Save per-image processed CSV
        out_csv = analysis_dir / (csv_path.stem + "_processed.csv")
        df.sort_values("rank").to_csv(out_csv, index=False)

        stem = csv_path.stem.replace("_oj_distribution", "")
        zrec = z_info.get(stem, {})
        summary_rows.append(
            {
                "File_Name": csv_path.name,
                "Number_of_Z_Stacks": zrec.get("number_of_z_stacks", "N/A"),
                "Z_Stack_Type": zrec.get("z_stack_type", "N/A"),
                f"Percentage_Aligned_within_{angle}°": pct_aligned,
                "Orientation_Mode": mode,
            }
        )
        logging.info(f"Processed CSV → {out_csv}")

    pd.DataFrame(summary_rows).to_csv(analysis_dir / "Alignment_Summary.csv", index=False)
    logging.info("Summary table written.")


# ────────────────────────────── Folder driver ────────────────────────────────

def process_one_folder(
    ij: "imagej.ImageJ",
    folder: str,
    fn_channel: int,
    angle: float,
    width: int,
    height: int,
) -> None:
    """Run Parts 1–3 for a single data folder."""
    angle_str = str(angle).replace(".", "_")
    res_dir, tbl_dir, img_dir = make_result_dirs(folder, angle_str)

    z_record = part1_z_projection(
        ij, folder, res_dir, fn_channel, width, height
    )
    part2_orientationj(res_dir, img_dir)

    analysis_dir = Path(res_dir) / "Analysis"
    analysis_dir.mkdir(exist_ok=True)
    part3_summarise(res_dir, str(analysis_dir), angle, z_record)


# ─────────────────────────────────── Main ────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description="Fibronectin alignment assay")
    ap.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input_paths.json that lists folders to analyse",
    )
    ap.add_argument("-a", "--angle", type=float, default=15.0, help="Angle threshold (°)")
    ap.add_argument("-w", "--width", type=int, default=500, help="Resize width (px)")
    ap.add_argument("-H", "--height", type=int, default=500, help="Resize height (px)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    ij = init_imagej()
    folders = read_folder_list(args.input)

    # channel index (1-based, as in ImageJ GUI)
    fn_ch = int(input("Enter fibronectin channel index (starting from 1): ").strip())

    ok = input("\nStart processing? (y/n): ").strip().lower()
    if ok not in ("y", "yes"):
        print("Cancelled by user."); return

    for fld in folders:
        process_one_folder(ij, fld, fn_ch, args.angle, args.width, args.height)

    print("\nAll folders processed. Disposing ImageJ…")
    ij.context().dispose()
    print("Done.")


if __name__ == "__main__":
    main()
