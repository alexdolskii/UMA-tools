# Active UMA-tools V2 workflow

Development for the forthcoming updated protocol is focused on this directory.
Other approaches in the repository are paused. V2 names the workflow under
development; the current Python package version is **0.2.6**.

See the [main README](../README.md) for installation, input JSON, parameters,
plate-template preparation, outputs, and updates.

## Commands and launchers

Run these stages individually in order with the same JSON. Before reporting,
place the plate template in the selected completed `Combined_Results` directory.

| Stage | Installed command | Launcher |
|---|---|---|
| Alignment | `uma_alignment` | `1_alignment.py` |
| Thickness | `uma_thickness` | `2_thickness.py` |
| Fibronectin area | `area_analysis` | `3_area.py` |
| Collect results | `uma_collect_results` | `4_collect_results.py` |
| Report | `uma_report` | `5_report.py` |

The numbered files are thin launchers for the same command implementations.
They accept identical arguments, for example, from the repository root:

```bash
python code/1_alignment.py -i input_paths.json -a 15
```

The old `alignment_analysis.py`, `thickness_analysis.py`, `area_analysis.py`,
`collect_results.py`, and `report.py` paths remain compatibility launchers.
Existing imports under `uma_tools` remain available. All launchers and installed
commands support `--help` and `--version` without starting Fiji.

## Module responsibilities

Reusable code belongs in the importable `uma_tools` package:

| Location | Responsibility |
|---|---|
| `uma_tools/cli.py` | Command arguments and completion status |
| `uma_tools/assays/` | Alignment, thickness, and area processing and measurements |
| `uma_tools/runtime/` | Shared headless Fiji initialization and worker cleanup |
| `uma_tools/collection.py` | Result selection and image correspondence |
| `uma_tools/reporting/` | Input validation, plate mapping, plots, and Excel output |
| `uma_tools/common/` | Configuration, file operations, output directories, logs, errors, version, and shared data definitions |

Reporting separates workflow, validation, plotting, and workbook generation.
Scientific dependencies load when needed. New code should import these focused
modules; compatibility launchers should contain no processing logic.

## Behavior to preserve

- Keep scientific calculations, operation order, thresholds, units, CSV
  columns, and report contents stable during structural changes. Numerical
  regression checks are required for processing changes.
- Use JSON `folder_paths` with absolute image-folder paths. For existing
  relative paths, area resolves against the JSON directory; the other stages
  resolve against the working directory. Shared helpers preserve this policy.
- Collection chooses the latest valid result for each assay independently.
  Reporting selects the latest completed collection, then requires all three
  summaries and one plate template there. Neither stage silently drops
  unmatched images; reporting does not fall back to an older collection.
- Preserve filename identities and explicit CSV/JSON writing policies.
  Exclude hidden files, macOS `._` files, and directories named like images.
- Give every run a new output directory and logs. Keep diagnostics with the
  run and close its log handlers before the next source folder.
- Dispose Fiji and close workers at the command boundary. Reusable processing
  functions must not terminate Python or close resources needed by later folders.

## Development checks

Run from the repository root after installing UMA-tools. Runtime requirements
are in [`environment_uma.yaml`](../environment_uma.yaml) and
[`pyproject.toml`](../pyproject.toml); the development checkers are separate:

```bash
python -m pip install -r requirements-dev.txt
ruff format --check code tests
ruff check code tests
pycodestyle --max-line-length=79 --max-doc-length=72 --ignore=E203,W503 code
python -m unittest discover -s tests -v
```

Use four-space indentation, `snake_case`, `UPPER_CASE` constants, 79-character
code lines, 72-character comments/docstrings, and explicit module interfaces.
The two pycodestyle exceptions match formatter slice spacing and PEP 8's
preferred break before binary operators. Required headless Matplotlib import
ordering uses a narrowly documented exception.

Enable the additional Fiji checks with:

```bash
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s tests -v
```

The [CI workflow](../.github/workflows/core-assays.yml) runs style, command,
numerical, and process-exit checks on Linux and macOS. Fiji tests use synthetic
TIFF data and require Java, Maven, and the initial component download.
Representative experimental ND2/TIFF data still need validation for their
channel order and calibration.
