# UMA-tools (Unit Matrix Assay)

Please visit full description of UMA-tools: [step-by-step-protocol](https://www.protocols.io/view/fibroblast-ecm-functional-units-a-medium-throughpu-gzpabx5if)

UMA-tools streamlines confocal analysis of 3D fibroblast/ECM units—turning time-consuming, error-prone steps into reproducible, scriptable workflows with clear, publication-ready outputs.

The project was developed in the [Edna (Eti) Cukierman lab](https://www.foxchase.org/edna-cukierman).

See also [fia-tools](https://github.com/alexdolskii/FIA-tools) for optimization of the fibroblast/ECM unit workflow.

For a complete guide to script usage, visit protocols.io.

# Usage
Automated image-analysis utilities for 3D fibroblast/ECM units assays in confocal microscopy.
For detailed description of 3D fibroblast/ECM units please refer [1](https://pubmed.ncbi.nlm.nih.gov/32222216/) and [2](https://pubmed.ncbi.nlm.nih.gov/27245425/).

In fibroblast/ECM 3D units, three readouts are especially informative: **fibronectin layer thickness**, **fibronectin fiber alignment**, and **nuclei counts/layering**. Together, these features function as rigorous quality-control metrics for unit maturation and as sensitive endpoints to quantify how compounds or genetic perturbations reshape the unit’s physiology (matrix organization, mechanics, and cellular architecture).

**Input** is a directory of .nd2 or .tif/.tiff confocal images representing Z-stacks with multiple detection channels (DAPI, fibronectin). To ensure correct channel mapping, keep the same channel order across every image.

## 1. Fibronectin Fiber Alignment — OrientationJ ImageJ/FIJI plugin (original protocol, Windows only)
The original script, Conda environment, and Windows setup instructions are collected in [Original Fibronectin Alignment Analysis](original_fibronectin_alignment_analysis/README.md).

The original fibronectin alignment protocol using the OrientationJ plugin was developed from the methods described in [this study](https://pubmed.ncbi.nlm.nih.gov/32222216/) and due to plugin constraints can be run only in ImageJ for Windows in GUI mode (headless mode is not supported).
- Pipeline: channel selection → 2D Max-Intensity Projection (XY) → resize/standardize → OrientationJ direction & coherence maps → orientation histograms and metrics.
- Outputs: color-coded orientation images, per-sample histograms, summary tables/Excel; % fibers within user-defined angle mode range.
- Windows

## 2. Fibronectin Fiber Alignment — orientationpy library (cross-platform)
Drop-in alternative to OrientationJ using a Python structure-tensor implementation (OrientationPy) for cross-platform, headless, and faster runs. Numeric values may differ slightly, but group-level trends remain consistent.
- Pipeline: projection & standardization → structure-tensor orientation + coherence → HSV orientation maps → orientation histograms/CSV.
- Outputs: identical report structure to the original protocol.
- Platforms: Linux, macOS.

## 3. 3D Unit Thickness Assay (fibronectin)
Using standard ImageJ/FIJI operations, we treat the fibronectin layer as a true 3D object (X–Y–Z): extract the FN channel, reslice to XZ to “look through” the layer, generate a 2D cross-sectional projection, enhance and background-correct, threshold to a binary mask, and compute Local Thickness—reporting the median thickness (with full distribution stats as optional outputs).
- Pipeline: channel selection → XZ reslice from 3D stacks → Max-Intensity Z-Projection → denoise (max filter + Gaussian) → background subtraction → Otsu threshold → Local Thickness (ImageJ plugin) → stats export.
- Outputs: per-image TIFF masks and thickness maps; CSV with area, mean/SD, min/median/max local thickness.
- Notes: assumes a consistent channel order across all images in a run.

## 4. Fibronectin Area — native-resolution SUM projection
The area assay projects the selected channel from original ND2/TIFF stacks into a 32-bit SUM image and applies an inclusive raw-intensity threshold. It reports FN-positive pixels, FN area, and the percentage of the full XY image occupied by FN. Original XY resolution and spatial calibration are retained. Each input folder is analyzed directly; filenames do not require a `_Seq####` identifier.

## 5. Nuclei Counts & Layer Prediction (3D)
The three analysis steps and their configuration are located in [alternative_nuclei_layers_assay](alternative_assays/alternative_nuclei_layers_assay). Each script uses the `nuclei_layers.json` next to it by default; `-i` selects a different configuration file.

Because fibroblast/ECM 3D units often exhibit strong background and debris that confound classical ImageJ thresholding, we perform StarDist 3D segmentation to robustly detect nuclei and extract their XYZ coordinates. We then apply scikit-learn spatial clustering to approximate nuclear “layers,” reporting both the layer count and per-nucleus membership. Note: for new cell types or staining conditions, you will likely need to train a custom StarDist model and point the script to it; step-by-step training and integration instructions are available on protocols.io.
- Pipeline: nuclei channel isolation → 3D denoising (Gaussian/mean) → StarDist 3D model (pre-trained models for fibroblastic lines; you may need to train your own) → QC overlays & tri-view projections → HDBSCAN clustering in 3D to infer layer-like groupings.
- Outputs: per-nucleus metrics (volume, centroid, equivalent diameter), image-level summaries, and study-level CSVs; QC figures for rapid validation.

## Alternative assays
Alternative nuclei-layer assays, marker-intensity analysis, data visualization tools, and StarDist models are grouped in [Alternative Assays](alternative_assays/README.md).

## Common features
- Headless batch processing across many folders/conditions listed in a single JSON file.
- User-guided channel selection (e.g., fibronectin, DAPI) with support for .nd2 and .tif/.tiff.
- Reproducible outputs: standardized images, per-image tables, and consolidated summaries suitable for downstream statistics and figure generation.

# Installation and commands for the main assays

The main programs in `code` use their own Conda environment, **uma_tools**, on macOS and Linux. Installing this repository registers `uma_alignment`, `uma_thickness`, `area_analysis`, `uma_collect_results`, and `uma_report` in that environment. The original Windows workflow and the tools in `alternative_assays` are separate approaches and are not included in this environment's supported scope.

## Install

Install [Git](https://git-scm.com/downloads) and [Miniforge](https://github.com/conda-forge/miniforge). On Apple Silicon, use the native arm64 installer; on an Intel Mac, use x86_64. Python and Java must use the same architecture.

Clone the development branch and run the installation commands from the repository root:

```bash
git clone --branch UMA-tools-V2 --single-branch https://github.com/alexdolskii/UMA-tools.git
cd UMA-tools
conda env create -f environment_uma.yaml
conda activate uma_tools
python -m pip install .
python -m pip check
```

For an existing UMA environment, update its dependencies and package using its actual name (for example, `uma_tools_new`), as described below. This does not require recreating the environment. To deliberately create a separate environment, use `conda env create -n uma_tools_v2 -f environment_uma.yaml` and activate that name before installing the package.

The environment specifies Python 3.10, OpenJDK 11, Maven, NumPy 1.26.4, PyImageJ 1.5.0, scyjava 1.10.0, jgo 1.0.4, and OrientationPy 0.3.0. scikit-image is installed automatically. Area uses the existing NumPy, PyImageJ, and scyjava dependencies. Reports add openpyxl 3.1.5 and Pillow >=10.1,<13 for Excel workbooks and embedded plots. TensorFlow and StarDist are not needed for these three assays. See [environment_uma.yaml](environment_uma.yaml) and [pyproject.toml](pyproject.toml) for the full requirements.

All three image-analysis programs initialize the fixed Fiji Maven endpoint `sc.fiji:fiji:2.14.0` in headless mode. The first analysis run downloads and caches Java components, including Fiji plugins, and requires access to Maven/SciJava repositories. A separate GUI installation of Fiji is not required for this route. The results collector uses only Python's standard library. Neither the collector nor the report command starts Fiji or asks for an image channel.

`uv.lock` records the Python dependency resolution. It does not install Java or Maven and does not replace the Conda setup above.

## Prepare the input JSON

Keep the UMA schema: the key is `folder_paths`, not FIA-tools' `paths_to_files`.

```json
{
  "folder_paths": [
    "/absolute/path/to/image folder"
  ]
}
```

Use absolute paths for image folders. The JSON file may be anywhere; pass its absolute path in quotes to run commands from any working directory. Hidden files, including macOS `._` metadata files, are excluded from image processing and file counts.

## Run

Activate the environment in each new terminal session:

```bash
conda activate uma_tools
uma_alignment --help
uma_thickness --help
area_analysis --help
uma_collect_results --help
uma_report --help
```

Fibronectin alignment, using the default 15-degree range:

```bash
uma_alignment -i "/absolute/path/input_paths.json"
```

To choose another angle, for example 10 degrees:

```bash
uma_alignment -i "/absolute/path/input_paths.json" -a 10
```

The program asks for the fibronectin channel (numbered from 1) and confirmation to start. It accepts `.nd2`, `.tif`, and `.tiff` images. Existing processing settings, including the 500-by-500 projection size and the 55% alignment classification threshold, are retained.

Fibronectin thickness:

```bash
uma_thickness -i "/absolute/path/input_paths.json"
```

The program asks for the file type (`.nd2` or `.tiff`), the fibronectin channel index, and confirmation. As in alignment, the channel prompt is `Enter fibronectin channel index (starting from 1):`. Enter `1` for a single-channel image. The channel count is read from each image, and the selected index is checked before extraction. These terminal questions are intentional; headless means no Fiji windows, not unattended execution.

Alignment and thickness retain their timestamped results directories and `log.log` files within the input folders. Alignment produces `Analysis/Alignment_Summary.csv`; thickness produces `Thickness_Summary.csv` with `Area`, `StdDev`, `Min`, `Max`, and `Median`. Image processing excludes macOS metadata files.

Fibronectin area, typically run after alignment and thickness:

```bash
area_analysis -i "/absolute/path/input_paths.json" -t 2000
```

The channel prompt is `Enter fibronectin channel index (starting from 1):`. It appears once for the entire command; enter `1` for a single-channel image. `--channel 3` supplies the index without a prompt. The actual channel count is checked in each original image.

The threshold syntax is `-t LOWER [UPPER]` or `--threshold LOWER [UPPER]`. Both endpoints are included in the mask:

| Arguments | Raw SUM intensity included in the mask |
|---|---|
| No `-t` | 2000 through the maximum finite float32 value |
| `-t 3000` | 3000 through the maximum finite float32 value |
| `-t 2000 50000` | 2000 through 50000 |
| `-t 2000 inf` | 2000 through the maximum finite float32 value |

The default upper bound is `3.4028234663852886e38`, not the maximum observed in one image. Bounds must be nonnegative, finite float32-representable values with UPPER >= LOWER; `inf` is accepted as an upper-bound shorthand. Effective float32 bounds and the original requested bounds are recorded. SUM also adds background and depends on Z count; the program does not normalize intensities or choose a threshold automatically.

Area reads every visible `.nd2`, `.tif`, and `.tiff` file directly in each JSON folder. It can run independently of alignment and thickness. Subfolders and previous result directories are excluded. `File_Name` and `Image_ID` retain the complete original filename, including its extension, so filenames without `_Seq####` and files with the same stem but different extensions remain distinguishable. ND2 files with multiple series, multiple time points, RGB-packed images, and invalid channels are rejected with diagnostics.

Every run creates `Area_assay_results_<timestamp>_<id>` inside each source folder. It contains `Projections_32bit`, `Masks`, `Fibronectin_Area_Summary.csv`, `FN_Projection_Manifest.csv`, `input_manifest.csv`, `run_parameters.json`, `run_status.json`, `run.log`, and `run_log.csv`. `FN_Area_Percent` uses the entire native-resolution XY frame as its denominator. `FN_Area` uses the original pixel calibration (square micrometers when calibrated in micrometers; square pixels for uncalibrated images). Masks are 8-bit TIFFs with values 0 and 255; saved SUM projections remain 32-bit float TIFFs.

An image error stops its source folder, retains partial CSVs plus `errors.csv` and `traceback.txt`, and allows the next JSON folder to run. Failures before processing create their own results directory in an available source folder, or in the working directory if no source folder is available. No logs are written beside the installed package. Area disposes ImageJ and closes the same worker pool as thickness before returning to the terminal. Any failed folder or shutdown error yields a nonzero process status.

The existing optional `--projection mean` / `--projection max` and `--projections-only` modes remain available. The default workflow is SUM with area measurement. `--projections-only` cannot be combined with `-t` and produces no masks or area summary.

### Collect the results of the three analyses

Use the same JSON after running alignment, thickness, and area:

```bash
uma_collect_results -i "/absolute/path/input_paths.json"
```

For **each source image folder separately**, the collector creates a new `Combined_Results_<source_folder_name>_<UTC_timestamp>` directory inside that folder. The name uses only the final folder component, never the full path. Repeated collections create separate directories. Unusual filename characters are replaced, and very long folder labels are shortened with a hash; the original folder name and full path remain in `run_status.json`.

The collector selects the latest valid result **independently for each analysis**, using the timestamp in the analysis directory name. It searches only the original analysis directories directly inside the JSON folder:

| Analysis | Input summary | Collected filename |
|---|---|---|
| Alignment | `Alignment_assay_results_angle_<angle>_<timestamp>/Analysis/Alignment_Summary.csv` | `<source_folder_name>_Alignment_Summary.csv` |
| Thickness | `Thickness_assay_results_<timestamp>/Thickness_Summary.csv` | `<source_folder_name>_Thickness_Summary.csv` |
| Area | `Area_assay_results_<timestamp>_<id>/Fibronectin_Area_Summary.csv` | `<source_folder_name>_Fibronectin_Area_Summary.csv` |

All collected CSVs retain their original bytes, column order, row order, and measurements. Tables are not joined. Older area results nested inside alignment directories, arbitrary renamed directories, symbolic-link result directories, hidden files/directories, `.partial.csv` files, and previous `Combined_Results` directories are excluded.

A valid result has a readable UTF-8 CSV, the assay's required columns, at least one image, unique image names, and finite values for its main measurements. Optional fields such as an unbounded area `Threshold_Upper` or unavailable alignment Z metadata may be blank or `N/A`. If `run_status.json` exists, it must report `SUCCESS`; recorded area image/mask counts must also agree with the summary row count. Newer invalid results are skipped with their reasons recorded, allowing an older valid result to be selected. Alignment and thickness do not have a machine-readable completion marker, so their final CSVs can be checked for structural validity but do not prove scientific correctness or coverage of every original image. Different valid runs with the same latest timestamp are reported as ambiguous instead of choosing arbitrarily.

**The main check is equality of the image sets across all selected, available tables.** Thickness and area retain complete original filenames. For alignment, only the exact terminal `_processed_orientation_distribution.csv` suffix is removed; the resulting stem is resolved against original filenames in the source folder and the selected tables. Filenames containing dots, spaces, Unicode, or no `_Seq####` identifier are supported. If `sample.nd2` and `sample.tif` could both correspond to one alignment row, the identity is ambiguous and collection fails. The check uses exact filenames, including case and extensions, not just row counts. It does not compare image contents, channels, or threshold settings; retain the selected assay CSVs and run directories to assess those settings.

If independently selected tables contain different images, the collector saves diagnostics only for that source folder. It does **not** search older runs for a combination that happens to match, remove unmatched rows, or create a successful collection of summary CSVs. Other JSON folders continue processing.

Missing analyses are allowed and reported in both the terminal and `run.log`, for example `Found 2 of 3 analyses`. The available tables must still match. With only one analysis, its CSV can be collected when its image identities are unambiguous, but the comparison is marked `NOT_COMPARABLE`. No empty CSV is created for a missing assay. No valid analyses, ambiguous identities, mismatches, or changed source CSVs produce a nonzero exit status.

Each collection starts with `run.log` and `run_status.json`. Selection adds `selection_report.csv` (selected source paths, timestamps, hashes, and reasons for skipping newer runs) and `image_check.csv` (per-image presence and correspondence). A successful run also contains the renamed summaries. Unexpected errors retain a traceback. Invalid configurations, including explicitly selected `._...json` files, are rejected before processing; if no source directory can receive diagnostics, a separate `Combined_Results` directory is created in the current working directory. Relative source paths, if used, are interpreted relative to the working directory, as in the existing assays; absolute paths are recommended.

The collector returns naturally to the terminal without starting Java. Exit status is `0` for successful collections, including collections with missing analyses, `1` if any source folder fails, and `130` if interrupted. Duplicate folder paths in the JSON are processed once.

### Generate the Excel report and plots

Place the experiment's annotation template (`.xlsx`) directly inside the completed `Combined_Results` directory, then run:

```bash
uma_report -i "/absolute/path/input_paths.json" --fn-threshold 20
```

Each JSON source folder is reported separately. The command selects the newest collector-completed `Combined_Results_<source_folder_name>_<timestamp>` directly inside that source folder, using the embedded timestamp. Eligible collector statuses are `SUCCESS` and `SUCCESS_WITH_MISSING_ANALYSES`; failed collector outputs are skipped. Ambiguous latest timestamps are errors. After selection, the chosen directory must contain all three collected assay CSVs and exactly one visible `.xlsx` annotation template directly inside it. A missing assay, missing or multiple templates, corrupt selected contents, or inconsistent image sets fails that source folder; the report does not fall back to an older collection or combine files from different collections. A collector status allowing missing analyses therefore does not waive the report's three-assay requirement.

The annotation template supplies the experiment's well and replicate information. Image correspondence uses complete original filenames, including extensions, and does not require `_Seq####`. Filenames must contain the supported `Well` token so images can be assigned to the template. Missing well annotations, duplicate or ambiguous image identities, and mismatches across assays produce diagnostics rather than dropping unmatched images. Hidden files, including macOS `._` files, are ignored.

`--fn-threshold` is a percentage from 0 through 100, with a default of 20. Images with `FN_Area_Percent` **strictly below** the threshold are excluded from the filtered report; images equal to the threshold are retained. Raw measurements remain available. This is a report filter on existing area percentages: the area assay's raw-intensity threshold, masks, and measurements are not recomputed.

Each successful run creates a new `UMA_Report_<source_folder_name>_<timestamp>` directory inside the selected `Combined_Results` directory. It includes:

- One `UMA_Report_<source_folder_name>_<timestamp>.xlsx` workbook with 20 sheets and 13 PNG plots, with the plots also embedded in the workbook.
- Raw, filtered, and excluded CSV tables and counts before and after filtering.
- Run logs, status, diagnostics, and a manifest recording the selected inputs and report settings.
- An `Inputs` directory preserving copies of the three source CSVs, annotation XLSX, `input_paths.json`, and `collector_run_status.json`, plus the collector's image check and selection report when available.

The report preserves the original plot semantics: points represent individual images, and colors distinguish technical replicates. Thickness `Area` is labeled in µm²; `StdDev`, `Min`, `Max`, and `Median` are labeled in µm. These labels assume the upstream thickness assay was run with micrometer calibration; report generation does not rescale measurements.

Failures retain diagnostics and allow subsequent JSON folders to run. Invalid JSON, including an explicitly selected `._...json`, is rejected before reading experiment data. If no source folder is available, failure logs are created in a separate report directory in the current working directory, never beside the installed package. The command returns naturally to the terminal without Java; a failed source folder produces a nonzero exit status.

The direct script interface remains available:

```bash
python code/alignment_analysis.py -i "/absolute/path/input_paths.json" -a 15
python code/thickness_analysis.py -i "/absolute/path/input_paths.json"
python code/area_analysis.py -i "/absolute/path/input_paths.json" -t 2000
python code/collect_results.py -i "/absolute/path/input_paths.json"
python code/report.py -i "/absolute/path/input_paths.json" --fn-threshold 20
```

To add the report command to an existing installation, update the `UMA-tools-V2` checkout and the existing environment's dependencies. Run these commands from the repository root and replace `uma_tools_new` with your actual environment name. An editable installation keeps subsequent Python source edits linked to the checkout:

```bash
git pull --ff-only
conda env update -n uma_tools_new -f environment_uma.yaml
conda activate uma_tools_new
python -m pip install --no-deps -e .
python -m pip check
uma_report --help
uma_report --version
```

The package version should be `0.2.5` or later; the report implementation version is `4.0.0`, the collector version is `1.0.0`, and the area version remains `2.1.0`. Registering a new command or changing package metadata still requires reinstalling the package, including for editable installs. Updating files with Git alone does not
replace a previously installed, non-editable package. Thickness runs report the
package version, implementation path, and reslice/projection calibration so the
installed implementation and spatial scale can be checked. For input calibrated
in micrometers, `Area` is in µm² and `StdDev`, `Min`, `Max`, and `Median` are in µm.
Rerun images processed with the uncalibrated projection; scaled filtering can
change the mask, so rescaling existing CSV values is not a general correction.

## Validation

Run the command/argument tests after installing the package:

```bash
python -m unittest discover -s tests -v
```

Results-collector tests cover latest-valid fallback, unavailable analyses, exact image identity checks, duplicate or ambiguous names, differing image sets with equal row counts, repeated and multi-folder collections, unchanged CSV bytes, and process exit without ImageJ. They use temporary files and require no experimental images or Java.

Report tests cover completed-collection selection, required inputs, exact image matching, plate annotations, threshold boundaries, corrected units, stable full/filtered plot positions and scales, workbook contents, archived input bytes, and error/interrupt handling. The installed-command test creates all 13 plots and verifies the 20-sheet workbook using synthetic measurements; it does not start Fiji.

Enable the additional Fiji integration tests explicitly:

```bash
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s tests -v
```

These integration tests process small synthetic TIFF stacks without a display, check output tables and metadata-file exclusion, compare calibrated masks, thickness maps, and measurements against the original projection macro, and compare direct Local Thickness results with the original menu command. Area command tests compare real SUM32 values and masks with known multichannel data, verify lower/upper threshold boundaries and calibrated area, check repeated runs without alignment folders, and exercise failure logging while continuing to another source folder. Separate child processes verify that thickness and area return to the terminal with the appropriate success/failure status. They require Java, Maven, and network access for the initial Fiji download. Synthetic TIFF tests do not establish accuracy on every experimental dataset; representative ND2 files and their channel/calibration metadata must also be validated on the target machine.

The [main-assay workflow](.github/workflows/core-assays.yml) builds the environment and runs these tests on Ubuntu and macOS. Check the latest GitHub Actions results before treating a platform as validated; a configured workflow alone is not evidence of a successful run.

## Other approaches

See [Original Fibronectin Alignment Analysis](original_fibronectin_alignment_analysis/README.md) for the Windows-only original protocol and its own environment, or [Alternative Assays](alternative_assays/README.md) for nuclei-layer analysis, marker-intensity analysis, visualization tools, and StarDist models.

# Dependencies and Tools Used
This program utilizes the following tools:

1. **Fiji** 
    This project used Fiji for preprocessing into preprocess image stacks as contrast enhancement, filtering, and particle analysis.

    [Fiji](https://fiji.sc/) is an open-source distribution of ImageJ focusing on image analysis. 
    
    - Repository: [Fiji](https://github.com/fiji/fiji)  
    - License: [GPL License](https://imagej.net/licensing/)

2. **StarDist**
    In this project, the standard StarDist model was employed to generate high-quality nuclei masks from image data, significantly improving segmentation accuracy and reducing background noise issues commonly encountered in immunofluorescence (IF) image analysis.
    
    [StarDist](https://stardist.net/)

    - Repository: [StarDist](https://github.com/stardist/stardist)  
    - License: [BSD 3-Clause License](https://github.com/stardist/stardist/blob/main/LICENSE.txt)
  
3. **OrientationPy**
   [OrientationPy](https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html) is a Python-based plugin used in this project for calculating the **alignment of fibronectin fibers**.

   - Repository: [OrientationPy](https://gitlab.com/epfl-center-for-imaging/orientationpy/)  
   - License: [The GNU General Public License, Version 3, 29 June 2007 (GPLv3)](https://gitlab.com/epfl-center-for-imaging/orientationpy/-/blob/main/LICENSE.md?ref_type=heads)

4. **OrientationJ plugin for ImageJ**:
   - Official page [OrientationJ](https://bigwww.epfl.ch/demo/orientation/)
   - Repository [Repository](https://github.com/Biomedical-Imaging-Group/OrientationJ)

5. **Skit-learn**
   - Official page [Skit-learn](https://scikit-learn.org/stable/)


## Contributors

- [Aleksandr Dolskii](mailto:aleksandr.dolskii@fccc.edu)

- [Ekaterina Shitik](mailto:shitik.ekaterina@gmail.com)

- Michael Miano


Enjoy your use 💫
