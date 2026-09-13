# Original Fibronectin Alignment Analysis

This folder contains the original UMA-tools approach for measuring fibronectin fiber alignment using Fiji/ImageJ and the OrientationJ plugin.

**Platform: Windows only.** This workflow requires an interactive Windows desktop session and runs Fiji/ImageJ in GUI mode. Headless execution is not supported.

## Files

- [alignment_analysis_original_approach.py](alignment_analysis_original_approach.py): the original fibronectin alignment analysis script.
- [environment_uma_original.yaml](environment_uma_original.yaml): the Conda environment definition for this approach, named `uma_original`.

## Analysis

The workflow generates 2D maximum-intensity projections, standardizes image size, and uses OrientationJ to calculate orientation maps and distributions. It measures the percentage of fibers within a selected angular window around the dominant direction, with a default window of +/-15 degrees. Images with an alignment percentage of at least 55% are classified as aligned.

Supported input formats are `.nd2`, `.tif`, `.tiff`, `.oif`, and `.oib`. Results include processed images, OrientationJ maps and distributions, normalized orientation images, and Excel/CSV analysis tables. Each run saves results in a timestamped `Alignment_assay_results_angle_<ANGLE>_<TIMESTAMP>` folder inside each input folder.

## Setup and launch on Windows

1. Place a Windows installation of Fiji with the OrientationJ plugin in `Fiji.app` at the repository root.
2. Edit [input_paths.json](../input_paths.json) at the repository root to list your image folders under `folder_paths`. Use forward slashes or escaped backslashes in JSON paths.
3. Open a Miniforge/Conda prompt in the repository root and create the environment:

```bat
conda env create -f original_fibronectin_alignment_analysis/environment_uma_original.yaml
conda activate uma_original
conda install -c conda-forge scikit-image
```

The script imports `scikit-image`, which is not explicitly listed in the original YAML; the additional installation command supplies this dependency. If `uma_original` already exists, skip the environment creation command and activate it.

Run the analysis from the repository root:

```bat
python original_fibronectin_alignment_analysis/alignment_analysis_original_approach.py
```

Follow the prompts to choose the angular window, specify the fibronectin channel, and confirm processing. The script reads `input_paths.json` and locates `Fiji.app` in the parent directory of this folder.
