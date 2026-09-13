# Alternative Assays

This directory groups alternative UMA-tools assays, analysis and visualization tools, reference data, and StarDist models.

| Folder | Contents |
| --- | --- |
| [alternative_nuclei_layers_assay](alternative_nuclei_layers_assay) | Three-step nuclei analysis: channel extraction, StarDist 3D segmentation, and quantification and layer clustering. Includes the JSON configuration and its loader. |
| [stardist_models_nuclei_layers](stardist_models_nuclei_layers) | Trained StarDist models, configuration and threshold files, model weights, and training logs. |
| [alternative_marker_intensity_analysis](alternative_marker_intensity_analysis) | R scripts for matrix-area quality control and pSMAD, PALLD, pFAK, and GS intensity analysis, with source tables and a conditions template in `DATA`. |
| [alternative_data_visualization](alternative_data_visualization) | File-renaming and unified plotting/statistical-analysis tools, JSON configurations, and usage notes. |

The original folder names and all files within them are preserved.

For nuclei analysis, edit [nuclei_layers.json](alternative_nuclei_layers_assay/nuclei_layers.json). Each of the three scripts loads this file from its own directory by default; use `-i` to select another configuration. For a bundled StarDist model, set `model_path` to the absolute path of the chosen directory containing `config.json` and the model weights under [stardist_models_nuclei_layers](stardist_models_nuclei_layers). Update any existing configuration that points to the models' former location.

See the [main README](../README.md) for commands run from the repository root.
