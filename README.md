# UMA-tools (Unit Matrix Assay)
UMA-tools streamlines confocal analysis of 3D fibroblast/ECM units—turning time-consuming, error-prone steps into reproducible, scriptable workflows with clear, publication-ready outputs.

The project was developed in the [Edna (Eti) Cukierman lab](https://www.foxchase.org/edna-cukierman).

For a complete guide to script usage, visit protocols.io.

# Aplication
Automated image-analysis utilities for 3D fibroblast/ECM units assays in confocal microscopy.
For detailed description of 3D fibroblast/ECM units please refer [1](https://pubmed.ncbi.nlm.nih.gov/32222216/) and [2](https://pubmed.ncbi.nlm.nih.gov/27245425/).

In fibroblast/ECM 3D units, three readouts are especially informative: **fibronectin layer thickness**, **fibronectin fiber alignment**, and **nuclei counts/layering**. Together, these features function as rigorous quality-control metrics for unit maturation and as sensitive endpoints to quantify how compounds or genetic perturbations reshape the unit’s physiology (matrix organization, mechanics, and cellular architecture).

**Input** is a directory of .nd2 or .tif/.tiff confocal images representing Z-stacks with multiple detection channels (DAPI, fibronectin). To ensure correct channel mapping, keep the same channel order across every image.

## 1. Fibronectin Fiber Alignment — OrientationJ ImageJ/FIJI plugin (original protocol, Windows only)
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

## 4. Nuclei Counts & Layer Prediction (3D)
Because fibroblast/ECM 3D units often exhibit strong background and debris that confound classical ImageJ thresholding, we perform StarDist 3D segmentation to robustly detect nuclei and extract their XYZ coordinates. We then apply scikit-learn spatial clustering to approximate nuclear “layers,” reporting both the layer count and per-nucleus membership. Note: for new cell types or staining conditions, you will likely need to train a custom StarDist model and point the script to it; step-by-step training and integration instructions are available on protocols.io.
Pipeline: nuclei channel isolation → 3D denoising (Gaussian/mean) → StarDist 3D model (pre-trained models for fibroblastic lines; you may need to train your own) → QC overlays & tri-view projections → HDBSCAN clustering in 3D to infer layer-like groupings.
Outputs: per-nucleus metrics (volume, centroid, equivalent diameter), image-level summaries, and study-level CSVs; QC figures for rapid validation.

## Common features
- Headless batch processing across many folders/conditions listed in a single JSON file.
- User-guided channel selection (e.g., fibronectin, DAPI) with support for .nd2 and .tif/.tiff.
- Reproducible outputs: standardized images, per-image tables, and consolidated summaries suitable for downstream statistics and figure generation.

# Installation 
For *2. Fibronectin Fiber Alignment — OrientationJ ImageJ/FIJI plugin* only:
- To download and install *OrientationJ plugin* please visit [OrientationJ](https://bigwww.epfl.ch/demo/orientation/)

For *all UMA-tools*:
- To download and install *git* please visit [Git Download page](https://git-scm.com/downloads).
- To download and install *conda* please visit [Miniforge github](https://github.com/conda-forge/miniforge)
- To download and install *OrientationPy* please visit [Official page: [Library description](https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html)] (version 0.3.0 is required)

## Usage
1. Before running the program, you need to modify a `input_paths.json` file. This file should contain a list of folders with .nd2 /.tiff/.tif images, and you can include as many folders as needed.
Additionally, before starting the program, make sure you know how many fluorescence channels you have (e.g., DAPI, Cy5) and their order in the file. You can check this by opening the image using the standard method in the GPU application (FiJi)[https://imagej.net/software/fiji/downloads].

2. Before first run only execute permission modofocation:
- Fibronectin Fiber Alignment — OrientationJ ImageJ/FIJI plugin (original protocol, Windows only)
```bash
chmod +x code/alignment_analysis_original_approach.py
```
- Fibronectin Fiber Alignment - orientationpy library (cross-platform)
```bash
   chmod +x code/alignment_analysis.py
```
-   3D Unit Thickness Assay (fibronectin)
```bash
  chmod +x code/thickness_analysis.py
```
-   Nuclei Counts & Layer Prediction (3D)
```bash
   chmod +x code/1_nla_fiji_channel_extraction.py
   chmod +x code/2_nla_stardist_prediction.py
   chmod +x code/3_nla_fiji_calculation.py
```

4. Run the main analysis script:
```bash
   - python code/alignment_analysis_original_approach.py -i input_paths.json
```
```bash
   - python code/alignment_analysis.py -i input_paths.json
```
```bash
   - python ccode/thickness_analysis.py -i input_paths.json
```
```bash
   - python code/1_nla_fiji_channel_extraction.py -i nuclei_layers.json
```
```bash
   - python code/2_nla_stardist_prediction.py -i nuclei_layers.json
```
```bash
   - python code/3_nla_fiji_calculation.py -i nuclei_layers.json
```



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
   - [Site]:(https://bigwww.epfl.ch/demo/orientation/)
   - [Repository]:(https://github.com/Biomedical-Imaging-Group/OrientationJ)

5. **Skit-learn**
   - [Site]:(https://scikit-learn.org/stable/)


## Contributors

- [Aleksandr Dolskii](aleksandr.dolskii@fccc.edu)

- [Ekaterina Shitik](mailto:shitik.ekaterina@gmail.com) 

Enjoy your use 💫

