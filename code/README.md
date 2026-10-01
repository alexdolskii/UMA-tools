# Active UMA-tools V2 workflow

Development for the forthcoming updated protocol is focused on this directory.
Other approaches in the repository are paused. V2 names the workflow under
development; the current Python package version is **0.2.11**.

See the [main README](../README.md) for installation, input JSON, parameters,
plate-template preparation, outputs, and updates.

## Five commands

Run the stages individually in order with the same JSON. Before reporting,
place the plate template in the selected completed
`<image_folder>/uma_assay/Combined_Results_.../` directory.
Keep original-image paths in the JSON. All five commands write results and
journals inside `uma_assay`; reports are saved beside collections. Collection
and reporting search only `uma_assay`, so rerun analyses saved in the earlier
layout.

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
Use the installed commands. Version 0.2.9 adds optional report statistics;
the image assay commands and measurements are unchanged. Runtime dependencies
remain unchanged in 0.2.11, which places core results and journals in `uma_assay`.

```bash
uma_report -i input_paths.json --fn-threshold 20 --stats-unit well
```

Omit `--stats-unit` to disable tests, or select `image` for exploratory tests
on individual images. Statistics use only FN-filtered data. Solid fills
define comparison blocks, and bold condition names identify each block's
control. See the [statistics rules](../README.md#optional-report-statistics)
for template validation, correction families, and interpretation.

## One package directory

The `code` directory contains this README and `uma_tools`, a package with
22 Python files. Its modules contain implementations rather than compatibility
adapters. `cli.py` routes commands directly to the five analysis/workflow
modules listed above; shared helpers are beside them in the same directory.

The remaining modules have these responsibilities:

| Module | Responsibility |
|---|---|
| `cli.py` | Command arguments, dispatch, and completion status |
| `__init__.py` | Package identity and installed version |
| `config.py` | JSON input folders and configuration errors |
| `files.py` | Core `uma_assay` location, filename labels, CSV/JSON writing, and checksums |
| `run.py` | Output directories, timestamps, and scoped logs |
| `progress.py` | Core command journals, archived logs, progress, and terminal prompts |
| `image_run.py` | Per-image outcomes and verifiable partial-run completion |
| `contracts.py` | Shared assay names, columns, and data definitions |
| `imagej.py` | Headless Fiji initialization and worker cleanup |
| `area_imagej.py` | SUM32 projection, threshold bounds, masks, and area measurements |
| `report_inputs.py` | Select, verify, and archive collected report inputs |
| `report_tables.py` | Read summary tables and the plate template |
| `report_validation.py` | Match images, validate annotations, and apply the FN coverage filter |
| `report_schema.py` | Report constants, data types, and validation errors |
| `report_statistics.py` | Validate color/bold controls and calculate optional Welch/Holm comparisons |
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
- Collection chooses the latest valid complete or audited partial result
  independently for each assay. Only registered failed images may be excluded
  from all collected copies; unexplained differences still fail validation.
  Reporting requires all three summaries and one plate template in the latest
  completed collection, without fallback. See the
  [partial-result rules](../README.md#progress-logs-and-partial-results).
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

The [test workflow](../.github/workflows/test.yml) runs numerical, report,
partial-result, logging and process-exit checks on Linux and macOS. The
[command workflow](../.github/workflows/core-assays.yml) checks installation
and lightweight entry points. Fiji tests use synthetic
TIFF data and require Java, Maven, and the initial component download.


Meaningful regressions cover latest-partial selection, unexplained mismatches,
unchanged source CSVs, excluded-image export, log rotation and routing,
invalid prompts, cancellation, NO_INPUT, and continuation after a damaged
image. Native checks also exercise all three installed image commands on a
valid TIFF and a damaged TIFF, verify Java logger restoration, and require
process exit within a timeout. Existing frozen numerical and calibration
checks remain in place.
