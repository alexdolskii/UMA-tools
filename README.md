# UMA-tools V2 (Unit Matrix Assay)

UMA-tools analyzes fibronectin alignment, layer thickness, and area coverage in
confocal images of 3D fibroblast/ECM units.

**Development status:** `UMA-tools-V2` is the second-version development branch
for a forthcoming updated protocol. Active development continues in
[`code`](code/README.md). Development of the other approaches is paused.
The current Python package version is **0.2.8**; this is separate from the V2
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
OrientationJ approach is linked below. Version 0.2.8 simplifies the package
into one directory with 19 Python files and removes duplicate launchers.

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
repository root directly into the latest completed `Combined_Results`
directory in each image folder. Fill in your experimental groups and save it.
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

All five commands support `--help` and `--version`. Run stages individually;
a command does not run earlier stages automatically. Module responsibilities
are documented in [`code/README.md`](code/README.md).

## Results and checks

Each run creates a separate timestamped directory with logs and diagnostics.
Source image folders are processed separately. Hidden files, including macOS
`._` files, are excluded from image processing and counts.

| Stage | Saved results |
|---|---|
| Alignment | Orientation images and tables; `Analysis/Alignment_Summary.csv` |
| Thickness | Masks, thickness maps, and `Thickness_Summary.csv` |
| Area | Native-resolution SUM32 projections, masks, and `Fibronectin_Area_Summary.csv` |
| Collection | `Combined_Results_<source_folder_name>_<timestamp>` inside each image folder; separate CSVs prefixed with that folder's name |
| Report | `UMA_Report_<source_folder_name>_<timestamp>` inside the selected collection; one Excel workbook, 13 plots, raw/filtered/excluded tables, and input copies |

Collection selects the newest valid result independently for each assay,
skipping newer invalid runs. **Selected tables must describe the same images**;
mismatches produce diagnostics instead of silently dropping rows. Missing
analyses are reported and may be collected, but reporting requires all three.
The report selects the newest completed collection and fails if its required
inputs are missing or inconsistent; it does not fall back to another collection.

For thickness calibrated in micrometers, `Area` is in µm² and `StdDev`, `Min`,
`Max`, and `Median` are in µm. Area coverage is a percentage of the full XY image.

## Update an existing installation

From your repository root on `UMA-tools-V2`, activate your existing environment
(replace `uma_tools_new` with its actual name):

```bash
git pull --ff-only
conda activate uma_tools_new
python -m pip install --no-deps -e .
python -m pip check
uma_alignment --version
```

An environment already set up for 0.2.5, 0.2.6, or 0.2.7 needs no dependency
changes for 0.2.8. The version command should report `uma_alignment 0.2.8`.
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
and [alternative assays](alternative_assays/README.md), including nuclei-layer
analysis, marker-intensity analysis, visualization tools, and StarDist models.
Their files remain available separately from the active V2 workflow.

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
