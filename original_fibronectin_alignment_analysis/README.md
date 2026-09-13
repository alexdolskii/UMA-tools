# Original Fibronectin Alignment Analysis

This folder contains the original UMA-tools approach for measuring fibronectin fiber alignment using Fiji/ImageJ and the OrientationJ plugin.

This folder can be used independently of the rest of the repository. The script reads `input_paths.json` and loads `Fiji.app` from its own directory, regardless of the current working directory.

**Platform: Windows only.** This workflow requires an interactive Windows desktop session and runs Fiji/ImageJ in GUI mode. Headless execution is not supported.

## Files

- [alignment_analysis_original_approach.py](alignment_analysis_original_approach.py): the original fibronectin alignment analysis script.
- [environment_uma_original.yaml](environment_uma_original.yaml): the Conda environment definition for this approach, named `uma_original`.
- [input_paths.json](input_paths.json): the local configuration listing the image folders to analyze.

## Analysis

The workflow generates 2D maximum-intensity projections, standardizes image size, and uses OrientationJ to calculate orientation maps and distributions. It measures the percentage of fibers within a selected angular window around the dominant direction, with a default window of +/-15 degrees. Images with an alignment percentage of at least 55% are classified as aligned.

Supported input formats are `.nd2`, `.tif`, `.tiff`, `.oif`, and `.oib`. Results include processed images, OrientationJ maps and distributions, normalized orientation images, and Excel/CSV analysis tables. Each run saves results in a timestamped `Alignment_assay_results_angle_<ANGLE>_<TIMESTAMP>` folder inside each input folder.

## Setup and launch on Windows

1. Install Fiji for Windows and the OrientationJ plugin. Place the complete `Fiji.app` directory inside `original_fibronectin_alignment_analysis`, next to the Python script. Fiji and OrientationJ are installed separately from the Conda environment.
2. Edit the [input_paths.json](input_paths.json) in this folder. Replace `/path/1` and `/path/2` with absolute paths to your image folders. Use forward slashes or escaped backslashes in JSON paths, for example:

```json
{
  "folder_paths": [
    "C:/Microscopy/experiment_1",
    "D:/Microscopy/experiment_2"
  ]
}
```

3. Open a Miniforge/Conda prompt and change to this folder, replacing the example path with its actual location:

```bat
cd /d "C:\path\to\original_fibronectin_alignment_analysis"
```

Create the environment on the first run:

```bat
conda env create -f environment_uma_original.yaml
conda activate uma_original
```

`scikit-image` is included in the YAML and is installed automatically with the `uma_original` environment.

If `uma_original` already exists, update it using the YAML in this folder instead of creating it again:

```bat
conda env update -n uma_original -f environment_uma_original.yaml
conda activate uma_original
```

Run the analysis from this folder:

```bat
python alignment_analysis_original_approach.py
```

Follow the prompts to choose the angular window, specify the fibronectin channel, and confirm processing. Keep the local `input_paths.json` and `Fiji.app` alongside the script when moving or copying this folder.
