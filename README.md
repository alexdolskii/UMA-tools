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

The main programs in `code` use their own Conda environment, **uma_tools**, on macOS and Linux. Installing this repository registers `uma_alignment`, `uma_thickness`, and `area_analysis` in that environment. The original Windows workflow and the tools in `alternative_assays` are separate approaches and are not included in this environment's supported scope.

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

For an existing UMA environment, activate its actual name (for example, `conda activate uma_tools_new`) and update the package as described below. Adding the area command does not require recreating the environment. To deliberately create a separate environment, use `conda env create -n uma_tools_v2 -f environment_uma.yaml` and activate that name before installing the package.

The environment specifies Python 3.10, OpenJDK 11, Maven, NumPy 1.26.4, PyImageJ 1.5.0, scyjava 1.10.0, jgo 1.0.4, and OrientationPy 0.3.0. scikit-image is installed automatically. Area uses the existing NumPy, PyImageJ, and scyjava dependencies. TensorFlow and StarDist are not needed for these three assays. See [environment_uma.yaml](environment_uma.yaml) and [pyproject.toml](pyproject.toml) for the full requirements.

All three programs initialize the fixed Fiji Maven endpoint `sc.fiji:fiji:2.14.0` in headless mode. The first analysis run downloads and caches Java components, including Fiji plugins, and requires access to Maven/SciJava repositories. A separate GUI installation of Fiji is not required for this route.

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

The direct script interface remains available:

```bash
python code/alignment_analysis.py -i "/absolute/path/input_paths.json" -a 15
python code/thickness_analysis.py -i "/absolute/path/input_paths.json"
python code/area_analysis.py -i "/absolute/path/input_paths.json" -t 2000
```

To install the area command and the previous thickness fixes, update the `UMA-tools-V2` checkout in your active UMA environment. An editable installation keeps subsequent Python source edits linked to the checkout:

```bash
git pull --ff-only &&
"$CONDA_PREFIX/bin/python" -m pip install --no-deps -e . &&
"$CONDA_PREFIX/bin/area_analysis" --version
```

The version should be `0.2.3` or later; the area script version is reported separately as `2.1.0`. Registering a new command or changing package metadata still requires reinstalling the package, including for editable installs. Updating files with Git alone does not
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
