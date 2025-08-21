# UMA-tools (Unit Matrix Assay)
UMA-tools streamlines confocal analysis of 3D fibroblast/ECM units—turning time-consuming, error-prone steps into reproducible, scriptable workflows with clear, publication-ready outputs.

The project was developed in the [Edna (Eti) Cukierman lab](https://www.foxchase.org/edna-cukierman).

For a complete guide to script usage, visit protocols.io.

## Aplication
Automated image-analysis utilities for 3D fibroblast/ECM (“unit”) assays in confocal microscopy. The toolbox focuses on fibronectin layer thickness, fibronectin fiber alignment, and 3D nuclei counting/layering. All workflows support batch, headless processing with ImageJ/FIJI and take a simple JSON manifest of input folders.

###
**3D Unit Thickness Assay (fibronectin)**
- Quantitatively measures the thickness of the fibronectin layer in 3D fibroblastic units.
- Pipeline: channel selection → XZ reslice from 3D stacks → Max-Intensity Z-Projection → denoise (max filter + Gaussian) → background subtraction → Otsu threshold → Local Thickness (ImageJ plugin) → stats export.
- Outputs: per-image TIFF masks and thickness maps; CSV with area, mean/SD, min/median/max local thickness.
- Notes: assumes a consistent channel order across all images in a run.

**Fibronectin Fiber Alignment — OrientationJ ImageJ/FIJI plugin (original protocol, Windows only)**
- Reproduces the original OrientationJ-based alignment analysis (OrientationJ is unreliable in headless mode cross-platform; use Windows).
Pipeline: channel selection → 2D Max-Intensity Projection (XY) → resize/standardize → OrientationJ direction & coherence maps → orientation histograms and metrics.
- Outputs: color-coded orientation images, per-sample histograms, summary tables/Excel; % fibers within user-defined angular windows (aligned vs disorganized categories).
- Platforms: Linux, macOS.

**Fibronectin Fiber Alignment — orientationpy library (cross-platform)**
Drop-in alternative to OrientationJ using a Python implementation (e.g., structure-tensor/gradient methods).
- Pipeline: projection & standardization → structure-tensor orientation + coherence → HSV orientation maps → orientation histograms/CSV.
- Outputs: identical report structure to the original protocol; small numeric differences may occur but group-level trends are consistent.
- Platforms: Linux, macOS.

**Nuclei Counts & Layer Prediction (3D)**
AI-assisted 3D nuclei segmentation and spatial clustering to approximate “layers”.
Pipeline: nuclei channel isolation → 3D denoising (Gaussian/mean) → StarDist 3D segmentation (pre-trained models for fibroblastic lines; you may need to train your own) → QC overlays & tri-view projections → HDBSCAN clustering in 3D to infer layer-like groupings.
Outputs: per-nucleus metrics (volume, centroid, equivalent diameter), image-level summaries, and study-level CSVs; QC figures for rapid validation.

**Common features**
- Headless batch processing across many folders/conditions listed in a single JSON file.
- User-guided channel selection (e.g., fibronectin, DAPI) with support for .nd2 and .tif/.tiff.
- Reproducible outputs: standardized images, per-image tables, and consolidated summaries suitable for downstream statistics and figure generation.

**Conventions**
One run can process multiple folders/conditions.
For reliable automation, keep fluorescence channel order identical across all images in a run.
If applying the nuclei pipeline to new cell types, plan to train a StarDist 3D model for best results.


## Installation 
To download and install *git* please visit [Git Download page](https://git-scm.com/downloads).

To download and install *conda* please visit [Miniforge github](https://github.com/conda-forge/miniforge)

Installing OrientationPy (v3+) The program requires a version higher than the one available via pip. Therefore, it needs to be installed via git (version 0.3.0 is used here). The official website is: https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html 

1. Clone the Repository. Download the  code for OrientationPy:
    - git clone https://gitlab.com/epfl-center-for-imaging/orientationpy.git
2. Enter the cloned directory:
    - cd orientationpy
3. Install OrientationPy and its dependencies:
    - pip install .
4. Return to the UMA-tools directory:
    - cd ..

## Usage
1. Before running the program, you need to modify a `input_paths.json` file. This file should contain a list of folders with .nd2 images, and you can include as many folders as needed.

Additionally, before starting the program, make sure you know how many fluorescence channels you have (e.g., DAPI, Cy5) and their order in the file. You can check this by opening the image using the standard method in the GPU application FiJi (https://imagej.net/software/fiji/downloads).


2. Run the main analysis script:
   - python ./code/alignment_analysis.py -i input_paths.json
   - python ./code/thickness_analysis.py -i input_paths.json



# Dependencies and Tools Used

This program utilizes the following tools:

1. **OrientationPy**  
   [OrientationPy](https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html) is a Python-based plugin used in this project for calculating the **alignment of fibronectin fibers**.

   - Repository: [OrientationPy](https://gitlab.com/epfl-center-for-imaging/orientationpy/)  
   - License: [The GNU General Public License, Version 3, 29 June 2007 (GPLv3)](https://gitlab.com/epfl-center-for-imaging/orientationpy/-/blob/main/LICENSE.md?ref_type=heads)

2. **Fiji**  
   [Fiji](https://fiji.sc/) (Fiji Is Just ImageJ) is an open-source distribution of ImageJ with a focus on image analysis. In this project, Fiji was used for **preprocessing image stacks**, including tasks such as contrast enhancement, filtering, and segmentation, to ensure the data is optimized for analysis in OrientationPy.

   - Repository: [Fiji](https://github.com/fiji/fiji)  
   - License: [GPL License](https://imagej.net/licensing/)


## References

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
   - Official page: [Library description](https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html)
   - Repository: [OrientationPy](https://gitlab.com/epfl-center-for-imaging/orientationpy/)


## Contributors

- [Aleksandr Dolskii](aleksandr.dolskii@fccc.edu)

- [Ekaterina Shitik](mailto:shitik.ekaterina@gmail.com) 

Enjoy your use 💫

