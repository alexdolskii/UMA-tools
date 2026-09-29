# Nuclei Layers Assay

Stage 1 is available as `uma_nla_prepare` in the existing UMA environment
on macOS and Linux. It uses UMA's pinned headless Fiji and closes its
background workers when the complete command finishes.

From the updated repository root, with the UMA environment active:

```bash
python -m pip install --no-deps ./nuclei_layers_assay
uma_nla_prepare --version
uma_nla_prepare -i nuclei_layers_assay/nuclei_layers.json
```

UMA-tools 0.2.9 must already be installed. This installs an additional
command; no environment recreation or new scientific dependencies are
needed. First complete the existing
[UMA Java/Fiji check](../README.md#verify-java-and-fiji-after-installation).

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

Stages 2 (`2_nla_stardist_prediction.py`) and 3
(`3_nla_fiji_calculation.py`) are unchanged and are **not yet integrated**
into the UMA environment. Stage 1 preserves their input directory and TIFF
naming conventions. Installing preparation does not install TensorFlow or
StarDist.

For development, install `pytest` separately and run
`python -m pytest nuclei_layers_assay/tests`. Real Fiji tests additionally
require `UMA_RUN_FIJI_TESTS=1`; they run in separate processes with a timeout
to check that the command actually returns to the terminal. To also check
a real single-field ND2 against the original filter sequence, set
`UMA_NLA_TEST_ND2=/absolute/path/to/file.nd2` (nuclei channel 1).
