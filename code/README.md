# Main UMA workflow

The numbered files document the order of the five commands. They are small
launchers; reusable code belongs in the importable `uma_tools` package.

| Stage | Launcher | Installed command |
|---|---|---|
| Alignment | `1_alignment.py` | `uma_alignment` |
| Thickness | `2_thickness.py` | `uma_thickness` |
| Fibronectin area | `3_area.py` | `area_analysis` |
| Collect results | `4_collect_results.py` | `uma_collect_results` |
| Report | `5_report.py` | `uma_report` |

Launchers accept the same arguments as the commands, for example:

```bash
python code/1_alignment.py -i input_paths.json -a 15
python code/2_thickness.py -i input_paths.json
python code/3_area.py -i input_paths.json -t 2000
python code/4_collect_results.py -i input_paths.json
python code/5_report.py -i input_paths.json --fn-threshold 20
```

The old `alignment_analysis.py`, `thickness_analysis.py`, `area_analysis.py`,
`collect_results.py`, and `report.py` filenames remain compatible launchers.
Existing imports under `uma_tools` are preserved through compatibility
modules. New internal code should import the focused modules below.

## Module responsibilities

| Location | Responsibility |
|---|---|
| `uma_tools/cli.py` | Parse command arguments and terminate each command with its result status |
| `uma_tools/assays/` | Image processing and scientific measurements for the three assays |
| `uma_tools/runtime/` | Shared headless Fiji initialization and background-thread cleanup |
| `uma_tools/collection.py` | Select completed assay results and compare image identities |
| `uma_tools/reporting/` | Validate collected tables and the plate map; generate plots and Excel |
| `uma_tools/common/config.py` | Reject metadata JSON files and apply explicit path-resolution policies |
| `uma_tools/common/files.py` | Write files with the requested CSV/JSON policies, hash files, and label folders |
| `uma_tools/common/run.py` | Allocate output directories and provide scoped and structured logs |
| `uma_tools/common/contracts.py` | Shared filename suffixes, columns, units, and result-directory grammar |
| `uma_tools/common/errors.py` | Validation failures with structured diagnostic details |
| `uma_tools/common/version.py` | Installed package version lookup |

The report has separate workflow, input validation, plotting, workbook,
dependency-loading, and data-contract modules. The renderer compatibility
facade exposes the existing API while these modules keep their individual
responsibilities. Scientific libraries load when needed; help and version
commands do not start Fiji.

## Compatibility rules

The layout change keeps command names, arguments, scientific calculations,
measurement units, thresholds, CSV columns, and report contents. Numerical
regression tests must pass before changing the processing modules.

Use one JSON with the `folder_paths` list and prefer absolute image-folder
paths. Relative-path policy remains explicit: area uses the JSON directory;
alignment, thickness, collection, and reporting use the working directory.
Shared configuration helpers must not silently change that distinction.

Other behaviors also remain deliberate:

- Collection chooses the newest valid result for each assay independently,
  including an older valid result when a newer result is invalid.
- Reporting chooses the newest successfully collected directory first. It
  then requires all three summaries and one plate template in that directory;
  it does not choose an older collection to compensate for missing inputs.
- Image identity checks retain full original filenames and reject ambiguous
  names. macOS `._` files, hidden files, and image-suffixed directories must
  not enter processing or image counts.
- Every run receives a separate results directory and appropriate logs.
  Failure diagnostics remain with that run. Scoped log handlers must close
  before processing the next folder.
- Fiji context disposal and worker-thread cleanup belong at the end of the
  command. Reusable processing functions must not terminate the Python
  process or shut down resources still needed by another folder.

CSV encodings, ordering, atomic writes, and error-status policies are explicit
parameters where the existing programs differ. Shared helpers must preserve
those contracts. External Java names and output table headers are public
interfaces; changing their spelling is not a style cleanup.

## Development checks

The runtime dependencies are defined in the root Conda environment and
`pyproject.toml`. The modular release does not change their versions. Install
the separate development requirements only when editing or checking code:

```bash
python -m pip install -r requirements-dev.txt
ruff format --check code tests
ruff check code tests
pycodestyle --max-line-length=79 --max-doc-length=72 --ignore=E203,W503 code
python -m unittest discover -s tests -v
```

Use four-space indentation, `snake_case` for Python functions and variables,
`UPPER_CASE` for constants, 79-character code lines, and 72-character comments
and docstrings. Keep function boundaries and inputs explicit, with type
annotations on shared contracts. Split functions by responsibility rather
than changing operation order to satisfy a size target.

Ruff checks imports, syntax, common errors, and formatting. Strict pycodestyle
also checks maintained code comments/docstrings. Its two exceptions match the
formatter's slice spacing (`E203`) and PEP 8's preferred break before binary
operators (`W503`). Required Matplotlib backend setup before `pyplot` imports
uses a narrowly documented exception; do not disable import-order checks
globally or eagerly load Java to satisfy the sorter.

The root GitHub Actions workflow runs formatting, style, launcher checks, and
the tests on Linux and macOS. Set `UMA_RUN_IMAGEJ_TESTS=1` to run the additional
Fiji integration tests locally; those tests require Java, Maven, and the
initial Fiji download. A passing style check does not replace numerical or
process-exit verification.
