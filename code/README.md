# Active UMA-tools V2 workflow

Development for the forthcoming updated protocol is focused on this directory.
Other approaches in the repository are paused. V2 names the workflow under
development; the current Python package version is **0.2.8**.

See the [main README](../README.md) for installation, input JSON, parameters,
plate-template preparation, outputs, and updates.

## Five commands

Run the stages individually in order with the same JSON. Before reporting,
place the plate template in the selected completed `Combined_Results` directory.

| Stage | Installed command | Implementation in `uma_tools` |
|---|---|---|
| Alignment | `uma_alignment` | `alignment_analysis.py` |
| Thickness | `uma_thickness` | `thickness_analysis.py` |
| Fibronectin area | `area_analysis` | `area_analysis.py` |
| Collect results | `uma_collect_results` | `collect_results.py` |
| Report | `uma_report` | `report.py` |

For example:

```bash
uma_alignment -i input_paths.json -a 15
```

All five commands support `--help` and `--version` without starting Fiji.
Version 0.2.8 removes duplicate script launchers and nested package layers.
Use the installed commands; their names, arguments, and calculations are
unchanged. Runtime dependencies are unchanged from versions 0.2.5–0.2.7.

## One package directory

The `code` directory contains this README and `uma_tools`, a package with
19 Python files. Its modules contain implementations rather than compatibility
adapters. `cli.py` routes commands directly to the five analysis/workflow
modules listed above; shared helpers are beside them in the same directory.

The remaining modules have these responsibilities:

| Module | Responsibility |
|---|---|
| `cli.py` | Command arguments, dispatch, and completion status |
| `__init__.py` | Package identity and installed version |
| `config.py` | JSON input folders and configuration errors |
| `files.py` | Filename labels, CSV/JSON writing, and checksums |
| `run.py` | Output directories, timestamps, and scoped logs |
| `contracts.py` | Shared assay names, columns, and data definitions |
| `imagej.py` | Headless Fiji initialization and worker cleanup |
| `area_imagej.py` | SUM32 projection, threshold bounds, masks, and area measurements |
| `report_inputs.py` | Select, verify, and archive collected report inputs |
| `report_tables.py` | Read summary tables and the plate template |
| `report_validation.py` | Match images, validate annotations, and apply the FN coverage filter |
| `report_schema.py` | Report constants, data types, and validation errors |
| `report_plots.py` | Generate the report figures |
| `report_workbook.py` | Build and verify the Excel workbook |

For custom Python integrations, import implementations directly, for example
`from uma_tools import alignment_analysis` or
`from uma_tools.report_plots import create_plots`. Earlier nested module paths
and compatibility adapters are no longer supported.

## Behavior to preserve

- Keep scientific calculations, operation order, thresholds, units, CSV
  columns, and report contents stable during structural changes.
- Use JSON `folder_paths` with absolute image-folder paths. Existing relative
  paths resolve against the JSON directory for area and against the working
  directory for the other stages.
- Collection chooses the latest valid result for each assay independently.
  Reporting selects the latest completed collection, then requires all three
  summaries and one plate template. Neither stage silently drops unmatched
  images; reporting does not fall back to an older collection.
- Preserve filename identities and explicit CSV/JSON writing policies.
  Exclude hidden files, macOS `._` files, and directories named like images.
- Give every run a new output directory and logs. Close its log handlers
  before the next source folder, and keep diagnostics with the run.
- Dispose Fiji and close workers at the command boundary. Reusable functions
  must not terminate Python or close resources needed by later folders.

## Development checks

Run from the repository root after installing UMA-tools. Runtime requirements
are in [`environment_uma.yaml`](../environment_uma.yaml) and
[`pyproject.toml`](../pyproject.toml); development checkers are separate:

```bash
python -m pip install -r requirements-dev.txt
ruff format --check code tests
ruff check code tests
pycodestyle --max-line-length=79 --max-doc-length=72 --ignore=E203,W503 code
python -m unittest discover -s tests -v
```

Use four-space indentation, `snake_case`, `UPPER_CASE` constants, 79-character
code lines, and 72-character comments/docstrings. The pycodestyle exceptions
match formatter slice spacing and PEP 8's break before binary operators.

Enable real Fiji checks with:

```bash
UMA_RUN_IMAGEJ_TESTS=1 python -m unittest discover -s tests -v
```

The [CI workflow](../.github/workflows/core-assays.yml) runs style, command,
numerical, and process-exit checks on Linux and macOS. Fiji tests use synthetic
TIFF data and require Java, Maven, and the initial component download.
