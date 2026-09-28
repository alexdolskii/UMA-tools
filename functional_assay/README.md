# Functional assay: stitching and cell analysis

Two commands for green-labelled cancer cells: stitch nine single-channel
ND2 Z-stacks per well, then measure the cell-mask area and count segmented
objects. Both commands use the same UMA environment and input JSON.

## Install and run

In your existing UMA environment, from the repository root:

```bash
conda activate uma_tools_new
python -m pip install --no-deps ./functional_assay
uma_stitching -i input_paths.json
uma_cell_count -i input_paths.json
```

Use your actual environment name. UMA-tools must already be installed in
that environment. No additional scientific dependencies are needed; follow
the main README's [Java/Fiji check](../README.md#verify-java-and-fiji-after-installation)
on a new computer. Both commands support `--help` and `--version` without
starting Fiji. To update, pull `UMA-tools-V2` and repeat the installation
command above; recreating the environment is unnecessary.

## 1. Stitching

The default overlap is **32.8%**. Override it with `--overlap 30`.
The JSON uses the existing `folder_paths` list:

```json
{"folder_paths": ["/absolute/path/to/original_images"]}
```

Each folder may contain several wells. Each well requires exactly nine
unique frames `0000`–`0008`, named like
`sample__WellB02_PointB02_0000_ChannelFITC_Seq0000.nd2`.
Each file must contain one channel. The fixed 3 × 3 grid is:

| Row | Left | Center | Right |
|---|---|---|---|
| Top | 0000 | 0001 | 0002 |
| Middle | 0005 | 0004 | 0003 |
| Bottom | 0006 | 0007 | 0008 |

The original Linear Blending settings and **Sharpen on all slices** are
retained. Output: `Stitched_Results/WellB02_stitched.tif` for each well.
**A new run deletes and replaces the previous `Stitched_Results`** when a
folder contains a valid well. Original ND2 files are retained.

Hidden/`._` files are ignored; explicitly selected `._…json` files are
rejected. Incomplete or duplicate frame sets are reported and skipped.
Fiji startup and worker shutdown use the existing UMA runtime. The Java
heap cap remains 16 GiB, as in the supplied script.

Each new run also writes `Stitched_Results/stitching_metadata.json`:
the overlap used, frame order and filenames, output dimensions and checksum,
physical pixel sizes from the nine original tiles, timestamp, and versions.
The stitching settings and image pixels are unchanged by this addition.

## 2. Cell count and mask area

`uma_cell_count` finds `Stitched_Results` inside each image folder in the
same JSON and reads only `Well…_stitched.tif` files. Keep the original nine
ND2 tiles in that folder: they are required to validate the physical scale.

Processing uses a MAX projection, rolling-ball background subtraction
(radius 50 pixels), Enhance Contrast (0.35% saturation, without intensity
normalization), and Smooth twice. The entire rectangular image is analyzed.

| Threshold option | Behavior |
|---|---|
| Omitted | Automatic `RenyiEntropy dark no-reset`, separately for each well |
| `-t` | Manual range 50–65535 |
| `-t 100` | Manual range 100–65535 |
| `-t 100 5000` | Manual range 100–5000 |

`--threshold` is equivalent to `-t`. Bounds apply to the processed 16-bit
projection and are inclusive; actual bounds are saved for every well.

The minimum particle area is **5 pixels²** by default. Set either
`--min-size-px 10` or `--min-size-um2 200`; these options are mutually
exclusive. The effective cutoff is recorded in both units.

- **Mask Area:** remove particles below the size cutoff, then measure the
  remaining positive pixels **before Watershed**. Large clusters remain
  included; holes are not filled. Edge-touching particles remain included.
- **Object Count:** apply Watershed to that cleaned mask, apply the size
  cutoff again, and exclude objects touching the rectangular image edge.
  The summed area of these counted objects is saved separately from Mask
  Area, because separation and edge exclusion can reduce it.

Areas are reported in **pixels² and µm²**. All nine ND2 files must contain
consistent, positive XY pixel sizes. A missing or inconsistent scale skips
that well with a logged error; other wells continue. Physical width and
height use the **actual stitched TIFF dimensions × original pixel size**.
Overlap is already represented in those dimensions and is not applied
again to area.

Every run creates a new folder alongside `Stitched_Results`:
`Cell_Analysis_<original-folder-name>_<UTC-timestamp>`. Previous cell-analysis
runs are retained. It contains:

- `Cell_Analysis_Summary.csv`: well, status, count, both area measures,
  thresholds, size cutoffs, physical scale, dimensions, and overlap.
- `Well…_objects.csv`: each counted object's area in both units.
- `Well…_area_mask.tif`: size-filtered mask before Watershed.
- `Well…_counting_mask.tif`: accepted objects after Watershed and filtering.
- `Well…_counted_contours.png`: green outlines on the processed MAX image.
- `run.log`, `run_status.json`, and a copy of `stitching_metadata.json`
  when available. TIFF masks contain physical XY calibration.

For older stitched TIFFs without metadata, the nine original ND2 files
still provide calibration; overlap is marked `not_recorded`, never assumed
to be 32.8%. New metadata is checked against the TIFF checksum before use.
To record overlap for an old dataset, rerun stitching with its intended
`--overlap` value (this replaces `Stitched_Results`).

A valid empty mask produces zero count/area and a header-only object CSV.
Failed wells have blank measurements and an error, rather than artificial
zeros. Exit code is `0` for a successful batch and `1` if any well/folder
fails; invalid command arguments use `2`. Fiji is disposed and its workers
are shut down before returning to the terminal.

## Checks

```bash
python -m unittest discover -s functional_assay/tests -v
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s functional_assay/tests -v
```

The second command also runs native Fiji comparisons and checks that a
Java worker process exits. Tests cover threshold modes, calibration from
nine tiles, size and edge filters, mask area before Watershed, physical
units, provenance, failure reporting, and preservation of earlier runs.
