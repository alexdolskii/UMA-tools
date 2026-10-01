# UMA-tools V2 (Unit Matrix Assay)

UMA-tools analyzes fibronectin alignment, layer thickness, and area coverage in
confocal images of 3D fibroblast/ECM units.

**Development status:** `UMA-tools-V2` is the second-version development branch
for a forthcoming updated protocol. Active development continues in
[`code`](code/README.md). Development of the other approaches is paused.
The current Python package version is **0.2.11**; this is separate from the V2
workflow name and the future protocol version.

## What changed

Compared with the [earlier workflow](https://github.com/alexdolskii/UMA-tools/tree/main):

| Aspect | Earlier workflow | Current V2 development |
|---|---|---|
| Launch | Individual scripts, such as `python code/alignment_analysis.py` | Installed commands in a dedicated `uma_tools` environment |
| Code | Processing organized mainly in individual scripts | Five commands backed by modules in one package directory |
| Main workflow | Alignment and thickness, alongside alternative assays | Alignment, thickness, FN area, result collection, and reporting |
| Results | Individual analysis outputs | Image-matched collected CSVs, plots, and an Excel report |

The modular code preserves the established calculations, units, parameters,
and output columns. Current alignment uses OrientationPy; the historical
OrientationJ approach is linked below. Version 0.2.9 adds optional comparisons
against plate controls and a filtered FN coverage plot. Image measurements
are unchanged.

## Install on macOS or Linux

Install [Git](https://git-scm.com/downloads) and
[Miniforge](https://github.com/conda-forge/miniforge). Use native arm64 packages
on Apple Silicon or x86_64 on Intel; Python and Java must share an architecture.

```bash
git clone --branch UMA-tools-V2 --single-branch https://github.com/alexdolskii/UMA-tools.git
cd UMA-tools
conda env create -f environment_uma.yaml
conda activate uma_tools
python -m pip install -e .
python -m pip check
```

To use another environment name, replace the creation command with
`conda env create -n uma_tools_new -f environment_uma.yaml` and activate
`uma_tools_new` instead. Requirements are defined in
[`environment_uma.yaml`](environment_uma.yaml) and [`pyproject.toml`](pyproject.toml).

The three image assays use Fiji `2.14.0` without graphical windows. The first
run downloads Java components and requires access to Maven/SciJava servers;
a separate Fiji GUI installation is unnecessary. Terminal questions still
appear. Collection and reporting do not start Fiji.

### Verify Java and Fiji after installation

**Complete this check before running image analyses on a new machine or after
updating the environment.** The tested UMA-tools setup uses **Fiji 2.14.0 with
OpenJDK 11**. Java and Maven are included in `environment_uma.yaml`; see the
[PyImageJ installation guidance](https://py.imagej.net/en/1.5.0/Install.html)
for background.

Activate the environment and check Java and Maven. Replace `uma_tools` with
your environment name if different, for example `uma_tools_new`.

```bash
conda activate uma_tools
java -version
mvn -version
```

Both commands should report **Java 11**. The patch version may differ, for
example `11.0.x`.

Next, verify that Python can initialize Fiji using the same startup procedure
as the analysis commands. Run this after installing the UMA-tools package
with `python -m pip install -e .` as shown above:

```bash
python - <<'PY'
from scyjava import jimport
from uma_tools.imagej import initialize_imagej, shutdown_imagej_workers

ij = None
try:
    ij = initialize_imagej()
    system = jimport("java.lang.System")
    print("Java version:", system.getProperty("java.version"))
    print("Java home:", system.getProperty("java.home"))

    if str(system.getProperty("java.specification.version")) != "11":
        raise RuntimeError("Use Java 11 for the tested UMA-tools environment.")
finally:
    try:
        if ij is not None:
            ij.dispose()
    finally:
        shutdown_imagej_workers()

print("Java/Fiji check passed.")
PY
```

The first initialization requires internet access to download Fiji components.
A successful check prints **`Java/Fiji check passed.`** and returns to the
terminal. This checks the Java runtime actually used by Python, in addition
to the Java executable found by the shell.

**If Java or Maven is missing**, install it in the activated UMA environment,
then reactivate that environment and repeat the checks:

```bash
conda install -c conda-forge "openjdk=11" maven
conda deactivate
conda activate uma_tools
```

Use your actual environment name in the last command.

**If these packages are installed but another Java version is reported**,
check environment activation, `JAVA_HOME`, and `PATH`. Python and Java must
also use compatible CPU architectures.

**If Fiji still fails to initialize**, run the
[PyImageJ diagnostic check](https://py.imagej.net/en/1.5.0/Troubleshooting.html):

```bash
python -c "import imagej.doctor; imagej.doctor.checkup()"
```

Include this output and the full initialization traceback, including preceding
Maven messages, when reporting the problem. The final
`Failed to create a JVM with the requested environment` message alone does
not identify the cause; dependency downloads, network access, or certificate
errors can also prevent initialization.

## Run the five stages

Use one JSON file with absolute paths to the original ND2/TIFF image folders:

```json
{
  "folder_paths": [
    "/absolute/path/to/image_folder"
  ]
}
```

Keep the same channel order across images. Activate the environment in each new
terminal session, then run the first four stages with the same JSON:

```bash
conda activate uma_tools
uma_alignment -i input_paths.json -a 15
uma_thickness -i input_paths.json
area_analysis -i input_paths.json -t 2000
uma_collect_results -i input_paths.json
```

Each image assay asks for the fibronectin channel, numbered from **1**.
Thickness also asks for the image file type. For execution outside the
repository, provide the JSON's absolute path in quotes.

Before the fifth stage, copy
[`UMA_96_well_plate_template.xlsx`](UMA_96_well_plate_template.xlsx) from the
repository root directly into the latest completed
`<image_folder>/uma_assay/Combined_Results_.../` directory.
Fill in your experimental groups and save it.
That directory must contain all three collected assay CSVs and **exactly one
plate-template `.xlsx`**. Then run:

```bash
uma_report -i input_paths.json --fn-threshold 20
```

The supplied template has one worksheet, `Plate Map`, and 96 empty well cells.
Keep columns 1–12 in `B1:M1` and rows A–H in `A2:A9`. Enter literal group names
in `B2:M9` (for example, well A01 is cell B2). Use identical spelling, case,
and spacing for wells in the same group. Do not use formulas or merge cells
in `A1:M9`. Every well represented in the images needs a group; unused wells
may stay blank. The empty template must be filled before reporting.
Image names must contain a supported well identifier such as `WellA02`;
`_Seq####` is not required.

**How the template is found:**

- The report selects the latest successful `Combined_Results` by the timestamp
  in its directory name, then searches directly inside that directory.
- The filename is unrestricted: `UMA_96_well_plate_template.xlsx`,
  `BK_far_day7.xlsx`, or a name with spaces all work. Renaming is optional;
  the name does not need to match the image folder or collected CSVs.
- Exactly one visible regular `.xlsx` file is required; `.XLSX` also works.
  Files starting with `.` (including `._`) or `~$`, symbolic links, directories,
  and Excel files inside subdirectories are ignored. `.xls` is not supported.
- No matching workbook causes an error. Two or more matching workbooks also
  cause an error, even if one is unrelated notes. The program does not choose
  by filename or switch to an older collection to find a template.
- The first worksheet in tab order is read by default. To choose another,
  use `--sheet "Worksheet name"` with its exact name; the active tab is not
  used to make this choice.

If a later collection run creates a newer successful `Combined_Results`, copy
the appropriate filled template into that new directory before reporting.

Key parameters:

- Alignment: `-a 15` sets the angle range to ±15 degrees, the default.
- Area: `-t 2000` sets the inclusive lower intensity threshold on the native
  32-bit SUM projection. The default upper bound is the largest finite float32
  value. Use, for example, `-t 2000 50000` to specify both bounds. Omitting `-t`
  uses 2000 as the lower bound.
- Report: `--fn-threshold 20` filters images with FN coverage **below 20%** from
  the filtered results; exactly 20% is retained. It does not recalculate masks.
- Report statistics are disabled unless `--stats-unit well` or
  `--stats-unit image` is supplied. The modes are described below.

All five commands support `--help` and `--version`. Run stages individually;
a command does not run earlier stages automatically. Module responsibilities
are documented in [`code/README.md`](code/README.md).

## Progress, logs, and partial results

Each command keeps a current log in
`<image_folder>/uma_assay/UMA_Logs/`:

| Step | Current log in `UMA_Logs/` |
|---|---|
| Alignment | `1_alignment.log` |
| Thickness | `2_thickness.log` |
| Area | `3_area.log` |
| Collection | `4_collect_results.log` |
| Report | `5_report.log` |

A new invocation moves that step's previous log into
`uma_assay/UMA_Logs/archive/`
with a timestamp. Other steps' logs stay in place. Detailed result-folder
logs are retained too. The terminal shows current activity and elapsed time,
with warnings and final counts; redirected output contains no animated control
codes. Full error traces belong in the logs. Prompts retry invalid answers;
`q`, EOF or Ctrl-C cancels with exit code 130. Help and version do not start
Fiji or create logs.

An image failure does not stop other images or source folders. Successful rows
remain available, while `run_status.json` records `PARTIAL`, every selected
filename, its outcome, and the failed stage/reason. `image_errors.csv` lists
failures. `SUCCESS` requires all selected images to finish; `FAILED` means
none succeeded, and `NO_INPUT` means no matching visible images were found.
Incomplete/error runs return a nonzero exit code, including partial runs.
A missing or unwritable source is reported and other folders continue. If no
source folder is available, startup diagnostics may use
`<working_directory>/uma_assay/`. They record the failure and do not create a
missing source folder.

**Collection selects the newest valid result for each assay, including an
audited `PARTIAL` run.** Its parameters take precedence over an older complete
run. Corrupt, cancelled, unfinished or invalid results are skipped with a
reason. A partial output without a verifiable image audit cannot justify
exclusions.

A registered failed image is excluded from **every collected CSV copy** so
that retained images match across analyses. Original analysis files are
preserved. Missing rows without a registered failure still fail validation;
the collector does not search for older runs just to obtain matching images.
An empty retained set fails collection. `processing_exclusions.csv`,
`image_check.csv` and `selection_report.csv` record the decisions and counts.
With no exclusions, collected summary bytes are unchanged.

Reporting still requires all three summaries and one plate template in the
latest successful collection. It uses the matched retained images, then
applies the requested FN filter. Processing failures are distinct from
low-FN exclusions: they appear in a separate **Processing Exclusions** Excel
sheet, a CSV and the report log. A collection missing an assay can be saved,
but cannot generate the three-assay report.

This change does not introduce process supervision or temporary-file cleanup.

## Results and checks

The five core commands write all analysis results, collected tables, reports
and journals inside `<image_folder>/uma_assay/`. The JSON continues to point
to the **original image folder**, not to `uma_assay`. Each source folder has
its own container, and each run gets a new timestamped directory.

Collection and reporting search **only this layout**, with no fallback to
results saved directly in the image folder. Rerun the image analyses and
collection when updating from the earlier layout. Existing results are not
moved or deleted. Image discovery reads only original files directly in the
source folder; it does not enter `uma_assay`. Hidden files, including macOS
`._` files, are excluded from image processing and counts.

| Stage | Saved results |
|---|---|
| Alignment | `uma_assay/Alignment_assay_results_angle_.../`: orientation images and tables; `Analysis/Alignment_Summary.csv` |
| Thickness | `uma_assay/Thickness_assay_results_.../`: masks, thickness maps, and `Thickness_Summary.csv` |
| Area | `uma_assay/Area_assay_results_.../`: native-resolution SUM32 projections, masks, and `Fibronectin_Area_Summary.csv` |
| Collection | `uma_assay/Combined_Results_<source_folder_name>_<timestamp>/`: separate CSVs prefixed with the original image folder's name; place one plate-template `.xlsx` here |
| Report | `uma_assay/UMA_Report_<source_folder_name>_<timestamp>/`: one Excel workbook, 14 plots, raw/filtered/excluded tables, input copies, and optional statistics |

Reports are siblings of collections. The report's `run_status.json` and input
provenance identify the selected `Combined_Results`; its template is still
read from that collection, regardless of the Excel filename.

Collection selects the newest valid result independently for each assay,
skipping newer invalid runs. **Selected tables must describe the same images**;
mismatches produce diagnostics instead of silently dropping rows. Missing
analyses are reported and may be collected, but reporting requires all three.
The report selects the newest completed collection and fails if its required
inputs are missing or inconsistent; it does not fall back to another collection.

For thickness calibrated in micrometers, `Area` is in µm² and `StdDev`, `Min`,
`Max`, and `Median` are in µm. Area coverage is a percentage of the full XY image.

## Optional report statistics

Without `--stats-unit`, reporting requires only the condition names in the
plate template. It produces seven full-data and seven FN-filtered plots,
without statistical tests, stars, or `ns`.

To enable comparisons, format the occupied well cells in the same template:

- Apply one identical **solid fill** to each comparison block, including its
  control. Exact stored color codes must match; theme colors and their tints
  are respected, so different shades form different blocks.
- Make the **entire condition name bold in every control well cell**. Each
  color must have exactly one control condition; its other conditions must
  not be bold. Cell positions do not define pairs.
- Use consistent color and bold formatting for all wells of one condition.
  Blank cells are ignored. Conditional formatting in the plate grid and
  partially formatted text are not supported. Invalid designs produce
  diagnostics instead of guessed controls.

```bash
uma_report -i input_paths.json --fn-threshold 20 --stats-unit well
```

`well` uses one mean of the retained images per well, with equal weight for
each well. Alternatively, `--stats-unit image` uses individual retained
images. This exploratory mode treats images as independent despite shared
wells; p-values can overstate evidence, and Holm does not correct that
dependence. Both modes describe technical comparisons within one plate,
not reproducibility across biological experiments.

Each noncontrol condition is compared only with its own control using a
two-sided Welch t-test. Tests use only FN-filtered data, for alignment, FN%,
and all five thickness metrics. Holm correction includes every planned
comparison across the seven metrics within each color block. For three
conditions versus one control, this is 21 hypotheses per block.

Fewer than two observations per arm in the selected mode, or undefined
variance, yields `Not tested` with a reason; unavailable tests remain in the
planned correction family. There is no automatic switch between units.
Filtered FN% results describe only images that passed the FN filter.

Filtered plots show adjusted significance: `*` for p < 0.05, `**` for p < 0.01,
`***` for p < 0.001, and `ns` otherwise. Points remain individual images;
captions identify the test unit and show image and well counts. Full-data
plots retain technical-well colors and red low-FN outlines. The selected
mode is recorded in the logs and run metadata.

Enabled statistics add `Well Means`, `Statistics`, and `Comparison Design`
Excel sheets and corresponding CSVs. The statistics include counts, means,
treatment-minus-control differences, nominal 95% Welch confidence intervals,
raw and adjusted p-values, and reasons for unavailable tests. Confidence
intervals are not adjusted for multiple comparisons. Differences in alignment
and FN% are expressed in percentage points.

## Update an existing installation

From your repository root on `UMA-tools-V2`, activate your existing environment
(replace `uma_tools_new` with its actual name):

```bash
git pull --ff-only
conda activate uma_tools_new
python -m pip install --no-deps -e .
python -m pip check
uma_report --version
```

An environment already set up for 0.2.5–0.2.10 needs no dependency changes for
0.2.11: no additional runtime dependencies are required. The version command
should report `uma_report 0.2.11`. Run the analyses again to populate the new
`uma_assay` layout before collecting results and generating reports.
For older environments, first update with
`conda env update -n uma_tools_new -f environment_uma.yaml`.
Recreating the environment is unnecessary. Editable installation (`-e`) makes
later Python source updates available after `git pull`; new commands or package
metadata changes still require reinstalling the package.

**Migration from earlier versions:** use the five installed commands above;
their names and arguments are unchanged. Duplicate script launchers have been
removed. Custom Python imports should use the module locations documented in
[`code/README.md`](code/README.md).

## Development and paused approaches

See [`code/README.md`](code/README.md) for module responsibilities and PEP 8
checks. The [automated workflow](.github/workflows/core-assays.yml) checks
commands, numerical regressions, and process completion on macOS and Linux.

Development is paused for the
[original Windows/OrientationJ approach](original_fibronectin_alignment_analysis/README.md)
and [alternative assays](alternative_assays/README.md), including
marker-intensity analysis and visualization tools.
Their files remain available separately from the active V2 workflow.

The separate [nuclei layers assay](nuclei_layers_assay/README.md) is being
integrated into the UMA environment: `uma_nla_prepare` and
`uma_nla_segment` are available. Its trained StarDist models now live in
the same folder. Stage 3 integration is still pending.

## Protocol and contributors

Developed in the [Edna (Eti) Cukierman lab](https://www.foxchase.org/edna-cukierman).
The [published protocol](https://www.protocols.io/view/fibroblast-ecm-functional-units-a-medium-throughpu-gzpabx5if)
describes the earlier workflow; the V2 update is in development for a future
protocol. Background methods: [reference 1](https://pubmed.ncbi.nlm.nih.gov/32222216/)
and [reference 2](https://pubmed.ncbi.nlm.nih.gov/27245425/).
Related project: [FIA-tools](https://github.com/alexdolskii/FIA-tools).

- [Aleksandr Dolskii](mailto:aleksandr.dolskii@fccc.edu)
- [Ekaterina Shitik](mailto:shitik.ekaterina@gmail.com)
- Michael Miano
