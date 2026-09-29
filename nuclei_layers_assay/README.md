# Nuclei Layers Assay

Stages 1 and 2 run in the existing UMA Python 3.10 environment on macOS
and Linux. Stage 1 (`uma_nla_prepare`) uses UMA's pinned headless Fiji and
closes its background workers when the command finishes. Stage 2
(`uma_nla_segment`) uses the local StarDist 3D models below; it does not
start Fiji.

## Install

From the updated repository root, activate the environment used for the
main UMA commands (replace `uma_tools_new` with your environment name):

```bash
conda activate uma_tools_new
python --version
python -m pip install "./nuclei_layers_assay[segment]"
uma_nla_prepare --version
uma_nla_segment --version
```

**Python must be 3.10.x**, as for UMA-tools. If it shows 3.12 or the prompt
shows `(base)`, activate the UMA environment before installing. UMA-tools
0.2.9 must already be installed.

The `[segment]` extra adds StarDist 0.9.1, CSBDeep 0.8.2 and TensorFlow 2.15.0.
It keeps UMA's NumPy 1.26.4 and requires no environment recreation. The five
bundled models are installed with the package, so the command also works
outside the repository. Model weights and thresholds are unchanged.
For preparation alone, `python -m pip install --no-deps ./nuclei_layers_assay`
remains available and does not install StarDist or TensorFlow.
First complete the existing
[UMA Java/Fiji check](../README.md#verify-java-and-fiji-after-installation).

## Stage 1: prepare

Edit the JSON to specify source folders and the nuclei channel (numbered
from 1). Use absolute folder paths; relative folder paths are interpreted
from the current working directory, as in the original NLA scripts.

```json
{
  "folder_paths": ["/absolute/path/to/images"],
  "nuclei_channel": 1,
  "gaussian_sigma": 4.0,
  "mean_radius": 3
}
```

The existing full `nuclei_layers.json` is also accepted. Settings for later
stages, such as `model_path`, are not needed for preparation.

```bash
uma_nla_prepare -i nuclei_layers_assay/nuclei_layers.json
```

The command processes regular ND2/TIF/TIFF files directly in each source
folder, excluding names beginning with `.` or `_`, including macOS `._`
metadata. A `._…json` configuration is rejected before reading. Each input
must contain a single field and a single time point; all Z slices of the
selected channel are retained.

Processing remains: **extract nuclei channel → 8-bit → Gaussian Blur 3D
→ Mean 3D**. Sigma and radius use pixel/voxel units, as in the original
script. Results keep their spatial calibration and are saved as
`processed/<original_stem>_nuclei.tif`. No projection or downsampling is
performed. `masks`, `analysis` and `clustering` are created for later stages.
Progress and per-image errors are recorded in `nuclei_analysis1.log` in
each source folder.

**Repeated runs:** as in the original script, preparation replaces all four
output folders: `processed`, `masks`, `analysis`, and `clustering`. Existing
results require one confirmation for the entire batch. For a deliberate
noninteractive replacement, add `--overwrite`. JSON, file-list and Fiji-startup checks
run before these folders are removed. An individual image failure
is logged; other images/folders are still processed and the command exits
with a nonzero status.

The previous direct launch remains available after installing this package:

```bash
python nuclei_layers_assay/1_nla_fiji_channel_extraction.py -i nuclei_layers_assay/nuclei_layers.json
```

## Stage 2: segment

Use the same JSON after preparation:

```bash
uma_nla_segment --list-models
uma_nla_segment -i nuclei_layers_assay/nuclei_layers.json
```

Set `model_path` to a short model name in the JSON, or override it with
`--model`. For example, **if this is the model trained for your data**:

```json
{
  "folder_paths": ["/absolute/path/to/images"],
  "nuclei_channel": 1,
  "gaussian_sigma": 4.0,
  "mean_radius": 3,
  "model_path": "stardist-512-v3-2",
  "n_tiles": [2, 4, 4]
}
```

The available names are `Stardist3D-512-60x`, `stardist-1024-v3`,
`stardist-256-v3`, `stardist-512-v3`, and `stardist-512-v3-2`. They now live
in [stardist_models_nuclei_layers](stardist_models_nuclei_layers). The nested
`model` directory of `Stardist3D-512-60x` is handled automatically. An
absolute custom model directory is also accepted. It must contain
`config.json`, `thresholds.json`, and `weights_best.h5` for a single-channel
3D model. Relative custom paths are checked against the working directory
and the JSON directory; ambiguous matches require an absolute path.

If `model_path` is `null` or absent, an interactive run asks once which
bundled model to use. Noninteractive runs require `model_path` or `--model`.
There is no automatic selection based on image dimensions or fallback to
`3D_demo`. Existing relative repository paths under the former model
location are mapped to the moved models when the old path no longer exists;
update old absolute paths or replace them with short names.

The command reads the regular TIFF stacks in each source folder's
`processed` directory. It excludes hidden/metadata files and rejects a
`._…json` configuration before reading. It retains the original global
1st/99.8th-percentile normalization, native XYZ dimensions and the model's
saved probability/NMS thresholds. `n_tiles` is ordered **Z, Y, X** and
controls inference tiling; it does not resize the data. The legacy
`downscale_factor` setting remains unused and is reported in the log.

Results are written to `masks`:

- `<prepared_stem>_mask.tif`: labelled instances, with the saved voxel scale.
  Labels use uint16 or uint32 when necessary to retain IDs above 65535.
- `<prepared_stem>_QC.png`: central XY/XZ/YZ slices and coloured instances.
- `segmentation_summary.csv`: per-image object counts and errors.
- `run_status.json`: completion status, model/thresholds, file checksums,
  library versions, tiling and per-image results.

`nuclei_analysis2.log` is written in each source folder. One model is loaded
for the entire batch. Individual image failures are logged and processing
continues; the command returns a nonzero exit status if any image fails.
Constant-intensity stacks receive an empty mask and a warning. Nonconstant
stacks with a zero percentile range are rejected for inspection.

**Repeated runs:** segmentation replaces `masks`, `analysis`, and
`clustering`, because previous calculations depend on the previous masks.
The `processed` input stacks are preserved. Existing nonempty results
require one confirmation for the batch; add `--overwrite` only for an
intentional replacement. Input/model checks and model loading run before
replacement. A failed image has no mask/QC output in the new run.

The previous launch form delegates to the same command after installation:

```bash
python nuclei_layers_assay/2_nla_stardist_prediction.py -i nuclei_layers_assay/nuclei_layers.json
```

Stage 3 (`3_nla_fiji_calculation.py`) is unchanged and is **not yet
integrated**. Stage 2 preserves its mask filenames and 3D label format.
Model compatibility and regression tests do not establish segmentation
accuracy for a new experiment; inspect the QC output for your chosen model.

## Development checks

Install `pytest` separately and run
`python -m pytest nuclei_layers_assay/tests`. Real Fiji tests additionally
require `UMA_RUN_FIJI_TESTS=1`; they run in separate processes with a timeout
to check that the command actually returns to the terminal. To also check
a real single-field ND2 against the original filter sequence, set
`UMA_NLA_TEST_ND2=/absolute/path/to/file.nd2` (nuclei channel 1).

Set `UMA_RUN_STARDIST_TESTS=1` for actual loading/prediction with all five
models and successful/failing process-exit checks. Set
`UMA_NLA_TEST_PREPARED=/absolute/path/to/prepared_nuclei.tif` as well for
exact mask comparison with the original stage-2 normalization and inference
on a representative prepared DAPI stack. The full comparison uses
`stardist-512-v3-2`; this is a regression fixture, not an automatic model
recommendation.
