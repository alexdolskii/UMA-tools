# Functional assay: stitching, cell measurements, and reporting

Stitch nine single-channel ND2 Z-stacks per well, then measure cell-mask
area and count segmented objects. Two reporting commands use these results:
`uma_functional_report` compares conditions at one time point;
`uma_survival_report` follows the same wells across explicit experiment days.
All four commands use the same UMA environment. The survival report uses
a separate JSON with the day assignments and one shared plate map.

## Install and run

In your existing UMA environment, from the repository root:

```bash
conda activate uma_tools_new
python -m pip install --no-deps .
python -m pip install --no-deps ./functional_assay
uma_stitching -i input_paths.json
uma_cell_count -i input_paths.json
# Add your Excel plate map to the selected Cell_Analysis_... folder first:
uma_functional_report -i input_paths.json --stats-unit well
```

Use your actual environment name. UMA-tools must already be installed in
that environment (UMA-tools 0.2.22 or later in the 0.2 series).
No additional scientific dependencies or environment recreation are needed; follow
the main README's [Java/Fiji check](../README.md#verify-java-and-fiji-after-installation)
on a new computer. All commands support `--help` and `--version` without
starting Fiji. To update, pull `UMA-tools-V2` and repeat the installation
commands above.

The reporting palette update is in `uma-functional-assay 0.5.1` and
`uma-tools 0.2.22`. Update **both** packages. Existing cell measurements
remain valid: rerun the reporting commands to regenerate the figures.

## Output locations, progress, and recovery

Keep JSON paths pointing to the **original image folders**. All new results
are stored inside `<original-images>/uma_functional_assay/`:

```text
uma_functional_assay/
  Stitched_Results/
  Cell_Analysis_<source>_<timestamp>/
  Functional_Report_<source>_<timestamp>/
  UMA_Logs/
    1_stitching.log
    2_cell_count.log
    3_functional_report.log
    archive/
```

Survival reports use `<output_dir>/uma_functional_assay/`, with
`Survival_Report_<experiment>_<timestamp>/` and `UMA_Logs/4_survival_report.log`.
Only this layout is searched. Results directly in original folders are
neither read nor moved; rerun stitching and cell analysis after this update.

The terminal shows the stage number and total, completed/total wells (or
days/plots), the current operation, and live elapsed time. A blocking Fiji
operation keeps its elapsed display active without inventing a percentage.
Redirected output uses plain progress lines. Detailed events, Python
warnings, and tracebacks are saved in each result's `run.log` and
`run_log.csv`, and in the numbered current journal. The previous journal is
archived when that command starts a new run for the source.

An unsuccessful well does not stop other wells or source folders.
`Processing_Exclusions.csv` records the well, stage, and reason; reports
also include a **Processing Exclusions** Excel sheet. `run_status.json`
records `SUCCESS`, `PARTIAL`, `FAILED`, `NO_INPUT`, or `CANCELLED`.
`PARTIAL` means usable results were saved with exclusions; it is not full
success. Unfinished `RUNNING` and cancelled runs are not report inputs.
Startup/output errors are recorded in the command journal when writable;
if no source folder is available, diagnostics use the working directory's
`uma_functional_assay/UMA_Logs/`.

Exit codes: **0** for complete success, **1** for partial/failed/no-input
results, **2** for invalid CLI arguments, and **130** for Ctrl-C. Reports
use the latest finalized `SUCCESS` or `PARTIAL` cell-analysis run. They
retain successful wells and never fill gaps from older runs. Missing values
stay missing; actual measured zeros remain valid.

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
retained. Output: `uma_functional_assay/Stitched_Results/WellB02_stitched.tif`
for each well. **A new attempt deletes and replaces the previous
`Stitched_Results`, including when the inputs are empty or invalid.** This
prevents a failed attempt from silently reusing stale TIFFs. Original ND2
files and timestamped cell-analysis/report folders are retained.

Hidden/`._` files are ignored; explicitly selected `._…json` files are
rejected. Incomplete or duplicate frame sets are reported and skipped.
Fiji startup and worker shutdown use the existing UMA runtime. The Java
heap cap remains 16 GiB, as in the supplied script.

Each new run also writes `Stitched_Results/stitching_metadata.json`:
the overlap used, frame order and filenames, output dimensions and checksum,
physical pixel sizes from the nine original tiles, timestamp, and versions.
The stitching settings and image pixels are unchanged by this addition.

## 2. Cell count and mask area

`uma_cell_count` finds `uma_functional_assay/Stitched_Results` inside each image folder in the
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
- `run.log`, `run_log.csv`, `run_status.json`, `Processing_Exclusions.csv`,
  and a copy of `stitching_metadata.json`. TIFF masks contain physical XY
  calibration. The summary is checkpointed after each well; its checksum
  binds the final CSV to the completion record.

Stitching metadata is required and checked against each TIFF checksum before
use. Only finalized `SUCCESS`/`PARTIAL` stitching runs are accepted. Wells
excluded by stitching stay excluded in cell analysis, with their reasons
preserved. Rerun stitching to use an older dataset in the new layout.

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

The program selects the **latest finalized `SUCCESS` or `PARTIAL`**
`Cell_Analysis_...` inside `<source>/uma_functional_assay/`. It uses the
timestamp in the folder name and verifies the counts and well outcomes in
`run_status.json`. Failed, unfinished, cancelled, or inconsistent completion
records are skipped with a message. A missing/invalid CSV or Excel file in the selected run
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

### Colors and well shapes in both reports

Boxplots use the same palette as `uma_report`, assigned separately within
each Excel color block. The Excel fill defines membership, not plot fill.
For up to **eight conditions, including the control**, the order is:

| Role | Color | HEX |
|---|---|---|
| Bold control | Grayish Lavender A | `#B5B1D8` |
| Treatment 1 | Dusky Green | `#004F46` |
| Treatment 2 | Orange | `#F37420` |
| Treatment 3 | Deep Indigo | `#051230` |
| Treatment 4 | Dull Blue Violet | `#80719E` |
| Treatment 5 | Ivory Buff | `#EBD3A2` |
| Treatment 6 | Violet | `#4F4086` |
| Treatment 7 | Verditter Blue | `#6FB5A8` |

Without an unambiguous bold control, descriptive plots use Dusky Green,
Grayish Lavender A, then the remaining colors in the table. Lavender then
represents an ordinary condition; it never creates a statistical control.
Statistics still require a valid marked control. For **nine or more
conditions**, treatments use Dusky Green tints; a marked control stays
lavender. Without a marked control, every condition uses a green tint.

Points are **neutral gray with dark outlines**. Shape identifies a technical
well within its condition: circle, square, triangle, then additional shapes
or numbers. Well order follows the full plate map (A01 through H12), so a
physical well keeps its shape and horizontal offset across endpoints,
days, baseline plots, and paired-change plots. Shapes restart in each
condition; the legend's Well 1, Well 2, etc. are positions within that
condition, not plate column numbers. Exact IDs are saved in `Well_Markers.csv`.

All annotated conditions and wells are assigned styles **before exclusions**.
Missing controls, wells, or complete days cannot change those assignments.
Repeated condition names in different Excel color blocks stay independent.
Survival panels add distinct hatches only with **six or more conditions**,
counting the full block, including control and conditions without data.
The first condition stays unhatched; the others use different patterns.
Single-time-point plots remain unhatched. Labels and shapes complement
colors; color alone is not a reliable identifier in grayscale.

Each report saves `Plot_Palette.csv`, `Well_Markers.csv`, and
`plot_palette.json`, plus **Plot Palette** and **Well Markers** Excel sheets.
The palette is a shared selection of historical digital colors, not a
numbered Wada combination; see the [main README](../README.md) for sources.
Measurements, masks, object contours, and statistical calculations are
unchanged by these display rules.

### Statistical comparisons

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
**beside the selected `Cell_Analysis_...`, inside `uma_functional_assay`**;
previous reports are retained. The Excel input stays inside `Cell_Analysis_...`.
It contains:

- One Excel workbook with the two plots, condition summaries, well data,
  all 96 well annotations/coverage, the plate map, and run details.
- `Object_Count.png` and `Mask_Area.png`.
- `Well_Data.csv`, `Condition_Summary.csv`, `Plate_Coverage.csv`, and
  `Processing_Exclusions.csv` (also in Excel).
- A `Statistics` Excel sheet and `Statistics.csv` only when requested.
- Exact input copies in `inputs/`, input checksums, plot provenance,
  `run.log`, and `run_status.json`.

The report validates well identifiers, completion status, numeric values,
and physical-area conversions, then reopens the workbook to verify saved
tables and plot embeddings. It needs **no Fiji, original ND2 files, or
stitched TIFFs**. Exit code is `0` when all reports succeed, `1` if any
folder fails, and `2` for invalid command arguments.

## 4. Survival assay over multiple days

Use `uma_survival_report` for **one plate imaged repeatedly**, with the same
well identifiers and conditions across days. Run stitching and cell analysis
with the existing `folder_paths` JSON first. No single-day report or local
Excel file inside each `Cell_Analysis_...` is required for this step.

Copy [survival_paths.example.json](survival_paths.example.json), rename it
`survival_paths.json`, and edit the paths and days:

```json
{
  "experiment_name": "survival_assay",
  "plate_template": "plate_map.xlsx",
  "output_dir": ".",
  "baseline_day": 1,
  "difference_days": [3, 5, 7],
  "timepoints": [
    {"day": 1, "folder": "day 1"},
    {"day": 3, "folder": "day 3"},
    {"day": 5, "folder": "day 5"},
    {"day": 7, "folder": "day 7"}
  ]
}
```

- Paths can be absolute or **relative to this JSON**, regardless of the
  terminal's working directory. In this example, place the JSON alongside
  the day folders and the Excel map. `output_dir: "."` saves reports in the
  `uma_functional_assay` subfolder there.
- `plate_template` is the exact path to your shared `.xlsx`, under any name.
  The grid, fill colors, and bold control rules are the same as in step 3.
  Use optional `"sheet": "Plate Map"` to select a worksheet; otherwise the
  first sheet is read. Instructions inside the older template about
  `Combined_Results` do not apply to this command.
- Each `folder` is an **original image folder**, not `Cell_Analysis_...`.
  Its latest finalized `SUCCESS`/`PARTIAL` analysis inside
  `uma_functional_assay` is selected independently. Missing/invalid data in
  that selected run make the day unavailable; older data are not substituted.
- Days are explicit, unique, non-negative integers, including **0**.
  Folder names do not need to contain day numbers. Plots are ordered by day.
- `baseline_day` must be a recorded day. `difference_days` explicitly lists
  which other days to subtract it from. The example produces `3-1`, `5-1`,
  and `7-1`; use baseline `0` and targets `[2, 3, 4]` for `2-0`, `3-0`, `4-0`.
  Reusing one folder for different days is rejected.

```bash
# Six descriptive plots and tables, without statistical testing:
uma_survival_report -i "/path/to/survival_paths.json"

# Add treatment-versus-control tests of the changes:
uma_survival_report -i "/path/to/survival_paths.json" --stats-unit well
```

For each **individual well**, the program calculates
`change = selected-day value - baseline-day value`. It retains negative
and zero changes. The two endpoints remain **Object Count** and **Mask Area
(µm²)** before Watershed. There is no percent-change or control normalization.

Missing measurements affect only the relevant day pair. A measured well
without a baseline stays in the raw table, but has no calculated change.
Each distribution shows its own `n`; missing measurements are never zeros.
If a later day is entirely missing or unusable, the report continues with
the other days and is marked **PARTIAL**. Its original measurements stay
absent, its changes stay blank, and untestable comparisons remain in the
planned Holm family. A baseline without usable wells mapped to the template
stops the report; day diagnostics are still saved. A partial baseline is
usable, with changes calculated only for wells present in both days.

**Additional/unannotated wells do not stop this multi-day report.** Their
measurements stay in `Raw_Measurements.csv` and Excel, marked `UNMAPPED`,
with warnings identifying the day and wells. They are excluded from grouped
plots and statistics because their condition is unknown. Annotated extra
wells can appear in raw distributions; changes still require both days.
The single-day `uma_functional_report` keeps its strict annotation check.

### Distributions and tests

Six PNGs are produced: Object Count and Mask Area, each shown **by day**,
for the **baseline alone**, and as **changes from baseline**. Each color
block has a panel, with condition distributions side by side within each
day, boxplots and all well points. Box colors identify conditions; hatches
add a second cue in blocks with six or more conditions. Point shapes identify
wells within each condition, consistently across all days and views. The
legend preserves full condition names and bold control. There are no line
plots. Raw and baseline plots are descriptive only.

With `--stats-unit well`, each treatment's **per-well changes** are compared
with the control's per-well changes using a two-sided Welch t-test.
**Holm correction covers both endpoints, all requested difference days,
and all planned control contrasts within each color.** For two treatments
and three day pairs, this is 12 tests per color. The same baseline can
contribute to several day pairs, but those days are never pooled as extra
replicates. Untestable comparisons remain in the planned family.

Only change plots carry `*`, `**`, `***`, or `ns`, based on adjusted p-values.
Fewer than two matched wells in either arm, or zero variance in both arms,
gives `Not tested`. The statistics table records the actual wells, counts,
mean changes, treatment-minus-control difference, unadjusted 95% confidence
interval, and raw/adjusted p-values. These describe technical variation
within this plate, not biological replication.

Baseline subtraction does not guarantee removal of differences in exposure,
focus, or segmentation. Processing parameters remain in the tables; changes
in scale, dimensions, overlap, size cutoffs, or threshold methods are logged.
Different automatic RenyiEntropy thresholds are retained as measured.

### Survival output

Each run creates `Survival_Report_<experiment_name>_<UTC-timestamp>` inside
`<output_dir>/uma_functional_assay`, preserving earlier reports. It contains:

- One Excel workbook with six embedded plots, raw measurements, per-well
  changes, condition summaries, coverage of all 96 wells for every day,
  selected analysis folders, the shared plate map, and run details.
- Six PNGs: `Object_Count_...` / `Mask_Area_...` with suffixes `By_Day`,
  `Baseline`, and `Changes`.
- `Raw_Measurements.csv`, `Changes_by_Well.csv`, `Group_Summary.csv`,
  `Well_Coverage.csv`, `Selected_Analyses.csv` (including each day's status
  and exclusion reason), and `Processing_Exclusions.csv` (day and well).
- `Statistics.csv` and the Excel Statistics sheet only when requested.
- Separate input copies for every day in `inputs/day_<number>/`, the JSON
  and plate map, SHA256 checksums, plot provenance, logs, and run status.

The report uses the saved summaries and completion records; it does not
start Fiji or reopen original images. Exit codes are `0` for success,
`1` for partial results or a report/input failure, `2` for invalid command
arguments, and `130` for interruption.

## Checks

```bash
python -m unittest discover -s functional_assay/tests -v
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s functional_assay/tests -v
```

The second command also runs native Fiji comparisons and checks that a
Java worker process exits. Tests cover threshold modes, calibration from
nine tiles, size and edge filters, mask area before Watershed, physical
units, provenance, failure reporting, and preservation of earlier runs.
Report tests also cover finalized partial-run selection, flexible Excel filenames,
missing annotations, independent color blocks, Welch/CI calculations and
Holm families, true zeros versus missing wells, optional statistics,
Excel round-trip checks, plotted well counts, and command startup without
Fiji. Synthetic test comparisons are not experimental results.
Survival tests additionally verify explicit day assignments (including day
zero), within-well subtraction, negative/zero changes, independent handling
of missing pairs, unmapped extra wells, Welch confidence intervals, Holm
across days, six complete plot exports, input preservation, and repeat runs.
Recovery checks cover a failed stitching well followed by successful wells,
propagation to cell analysis, stale-TIFF rejection, ignored old locations,
per-command journal rotation, live progress shutdown, folder I/O failures,
missing whole days, absent/unmapped baselines, audited exclusions, and
reduced sample counts without substituting zero or older measurements.
