"""
Compatibility API delegating to independent report components.

Importing this facade does not load NumPy, Matplotlib, Pillow, or
Openpyxl.
"""

from .constants import (
    ALIGNMENT_FILENAME_SUFFIX,
    ALIGNMENT_PATTERN,
    BASE_COLORS,
    EVENT_COLUMNS,
    FN_INCLUDED_FLAG,
    FN_LOW_FLAG,
    FN_METRIC,
    FN_REASON_COLUMN,
    FN_THRESHOLD_COLUMN,
    LOW_FN_EDGE_COLOR,
    NUMBER_PATTERN,
    POINT_PATTERN,
    POINT_WELL_PATTERN,
    SCRIPT_VERSION,
    SEQUENCE_PATTERN,
    SHEET_NAMES,
    THICKNESS_METRICS,
    THICKNESS_UNITS,
    WELL_PATTERN,
)
from .dependencies import (
    load_dependencies,
)
from .io import (
    RunLog,
    save_csv,
    save_details,
    save_json,
    sha256_file,
    utc_now,
)
from .models import (
    EventLogger,
    ParsedRecord,
    ReportData,
    ReportInputs,
    ValidationError,
)
from .plate import read_template
from .plots import (
    box_definition,
    create_plots,
    replicate_colors,
)
from .source_tables import (
    metadata_column,
    numeric_values,
    optional_number_sort,
    original_stem_lookup,
    parse_filenames,
    read_csv_table,
    validate_fibronectin,
    validate_fn_threshold,
)
from .validation import filter_counts, validate_and_merge
from .workbook import (
    build_workbook,
    put_cell,
    title_sheet,
    verify_workbook,
    write_table,
)

__all__ = [
    "SCRIPT_VERSION",
    "THICKNESS_METRICS",
    "THICKNESS_UNITS",
    "FN_METRIC",
    "FN_THRESHOLD_COLUMN",
    "FN_LOW_FLAG",
    "FN_INCLUDED_FLAG",
    "FN_REASON_COLUMN",
    "LOW_FN_EDGE_COLOR",
    "ALIGNMENT_PATTERN",
    "ALIGNMENT_FILENAME_SUFFIX",
    "SEQUENCE_PATTERN",
    "WELL_PATTERN",
    "POINT_WELL_PATTERN",
    "POINT_PATTERN",
    "NUMBER_PATTERN",
    "BASE_COLORS",
    "SHEET_NAMES",
    "EVENT_COLUMNS",
    "load_dependencies",
    "RunLog",
    "utc_now",
    "save_json",
    "save_csv",
    "save_details",
    "sha256_file",
    "ValidationError",
    "ReportData",
    "ReportInputs",
    "EventLogger",
    "ParsedRecord",
    "read_csv_table",
    "parse_filenames",
    "original_stem_lookup",
    "optional_number_sort",
    "read_template",
    "numeric_values",
    "metadata_column",
    "validate_fn_threshold",
    "validate_fibronectin",
    "filter_counts",
    "validate_and_merge",
    "replicate_colors",
    "box_definition",
    "create_plots",
    "put_cell",
    "write_table",
    "title_sheet",
    "build_workbook",
    "verify_workbook",
]
