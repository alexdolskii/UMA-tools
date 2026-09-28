# Functional assay: stitching, cell analysis, and plate reports

Three commands for green-labelled cancer cells: stitch nine single-channel
ND2 Z-stacks per well, measure cell-mask area and count segmented objects,
then compare conditions using a 96-well plate map. All three commands use
the same UMA environment and input JSON.

## Install and run

In your existing UMA environment, from the repository root:

```bash
conda activate uma_tools_new
python -m pip install --no-deps ./functional_assay
uma_stitching -i input_paths.json
uma_cell_count -i input_paths.json
# Add your Excel plate map to the completed Cell_Analysis_... folder first:
uma_functional_report -i input_paths.json --stats-unit well
```

Use your actual environment name. UMA-tools must already be installed in
that environment. No additional scientific dependencies are needed; follow
the main README's [Java/Fiji check](../README.md#verify-java-and-fiji-after-installation)
on a new computer. All commands support `--help` and `--version` without
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

## 3. Plate report and control comparisons

Each source folder in the JSON represents **one plate at one time point**.
`uma_functional_report` creates a separate report for each source folder;
it does not combine plates, time points, or repeated analysis runs.

```bash
# Tables and the two endpoint plots, without statistical tests:
uma_functional_report -i input_paths.json

# Add well-level comparisons against each color block's control:
uma_functional_report -i input_paths.json --stats-unit well
```

### Select the results and add the plate map

The program selects the **latest fully completed** `Cell_Analysis_...`
directly inside each source folder. It uses the timestamp in the folder
name and verifies `run_status.json`. Newer failed or partial runs are
skipped with a message. A missing CSV or Excel file in the selected run
is an error: the program does **not** fall back to an older marked run.

Place exactly **one `.xlsx` directly inside that selected folder**. Its
filename is arbitrary; renaming it does not change discovery. Hidden,
`._`, and Excel lock files (`~$...`) are ignored, as are nested folders.
The first worksheet is used by default; select another with
`--sheet "Plate Map"`.

Use the repository's [96-well template](../UMA_96_well_plate_template.xlsx):
keep rows **A-H** in `A2:A9`, columns **1-12** in `B1:M1`, and enter conditions
in `B2:M9`. For this functional assay, use the `Cell_Analysis_...` location
and command above, rather than the template's older `Combined_Results` /
`uma_report` instructions.

- **Cell text** identifies the condition, e.g. `Control` or `Treatment 1`.
- **Solid fill color** defines a comparison block. Use exactly the same
  color code within a block; different shades/tints are separate blocks.
- **Bold text on the whole cell** marks the control condition. With
  statistics enabled, each block must have exactly one control condition
  and at least one treatment. Mark all wells of that control consistently.
- A condition name such as `Control` may recur in different colors; those
  conditions remain separate. Within one color, identical text pools wells.
- Use literal names and direct formatting. Formulas, merged grid cells,
  conditional formatting, mixed bold within a cell, and inconsistent bold
  among wells of the same condition are rejected. Avoid surrounding spaces.

A measured well with a blank plate-map cell **stops the report** and is
identified by well and Excel coordinate in the console, log, and
`validation_errors.csv`. Annotated wells without measurements are listed
as missing; they never become zero values. Actual measured zeros are kept.
Without statistics, missing/ambiguous controls are reported as warnings;
the measurements can still be plotted.

### Measurements, plots, and statistics

There are exactly **two endpoint plots**: **Object Count** and **Mask Area
(µm²)**. Mask Area uses the size-filtered mask **before Watershed**, as in
step 2. Other measurements and processing parameters remain in the tables;
no FN filter or normalization to control is applied.

Each plot has a panel per color block, with all measured wells as points,
the number of wells (`n`), and a boxplot where at least two wells are
available. The box shows median and interquartile range; whiskers extend
to the most extreme observed values within 1.5 times that range. The
control label is bold. A condition without measured wells shows `n=0`.

Statistics run **only with `--stats-unit well`**. One well is one technical
replicate; the nine tiles, individual cells, and repeated processing are
not additional replicates. Each treatment is compared pairwise with its
block's control using a two-sided Welch t-test. Holm correction covers
**both endpoints and all planned control comparisons within that color**.
Separate colors are separate correction families.

Graphs show adjusted-p labels: `* <0.05`, `** <0.01`, `*** <0.001`, or
`ns >=0.05`. Fewer than two wells in either arm, or zero variance in both
arms, produces **Not tested**, never `ns`. Unavailable comparisons remain
in the planned Holm family. Excel and CSV record well counts, means,
treatment-minus-control differences, **unadjusted** 95% confidence
intervals, raw/adjusted p-values, and reasons for untested comparisons.
These are within-plate technical comparisons, not biological replication.

### Saved report

Each run creates `Functional_Report_<original-folder-name>_<UTC-timestamp>`
**inside the selected `Cell_Analysis_...`**; previous reports are retained.
It contains:

- One Excel workbook with the two plots, condition summaries, well data,
  all 96 well annotations/coverage, the plate map, and run details.
- `Object_Count.png` and `Mask_Area.png`.
- `Well_Data.csv`, `Condition_Summary.csv`, and `Plate_Coverage.csv`.
- A `Statistics` Excel sheet and `Statistics.csv` only when requested.
- Exact input copies in `inputs/`, input checksums, plot provenance,
  `run.log`, and `run_status.json`.

The report validates well identifiers, completion status, numeric values,
and physical-area conversions, then reopens the workbook to verify saved
tables and plot embeddings. It needs **no Fiji, original ND2 files, or
stitched TIFFs**. Exit code is `0` when all reports succeed, `1` if any
folder fails, and `2` for invalid command arguments.

## Checks

```bash
python -m unittest discover -s functional_assay/tests -v
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s functional_assay/tests -v
```

The second command also runs native Fiji comparisons and checks that a
Java worker process exits. Tests cover threshold modes, calibration from
nine tiles, size and edge filters, mask area before Watershed, physical
units, provenance, failure reporting, and preservation of earlier runs.
Report tests also cover completed-run selection, flexible Excel filenames,
missing annotations, independent color blocks, Welch/CI calculations and
Holm families, true zeros versus missing wells, optional statistics,
Excel round-trip checks, plotted well counts, and command startup without
Fiji. Synthetic test comparisons are not experimental results.
