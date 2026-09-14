"""Report measurements, output schema, and validation data contracts."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Protocol, TypedDict

from .contracts import ALIGNMENT_SUFFIX
from .contracts import EVENT_COLUMNS as EVENT_COLUMNS
from .contracts import THICKNESS_METRICS as THICKNESS_METRICS
from .contracts import THICKNESS_UNITS as THICKNESS_UNITS

SCRIPT_VERSION = "4.0.0"
FN_METRIC = "FN_Area_Percent"
FN_THRESHOLD_COLUMN = "FN_Area_Threshold_Percent"
FN_LOW_FLAG = "Below_FN_Threshold"
FN_INCLUDED_FLAG = "Included_In_Filtered_Plots"
FN_REASON_COLUMN = "Exclusion_Reason"
LOW_FN_EDGE_COLOR = "#D62728"
ALIGNMENT_PATTERN = re.compile(
    r"^Percentage_Fibers_Aligned_Within_([0-9]+(?:\.[0-9]+)?)_Degree$"
)
ALIGNMENT_FILENAME_SUFFIX = ALIGNMENT_SUFFIX
SEQUENCE_PATTERN = re.compile(r"(?:^|[_. -])Seq([0-9]+)(?=[_. -]|$)")
WELL_PATTERN = re.compile(
    r"(?:^|[_. -])Well([A-Ha-h])(0?[1-9]|1[0-2])(?=[_. -]|$)"
)
POINT_WELL_PATTERN = re.compile(
    r"(?:^|[_. -])Point([A-Za-z])([0-9]+)(?=[_. -]|$)"
)
POINT_PATTERN = re.compile(
    r"(?:^|[_. -])Point([A-Ha-h])(0?[1-9]|1[0-2])_([0-9]+)(?=[_. -]|$)"
)
NUMBER_PATTERN = re.compile(
    r"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"
)
BASE_COLORS = [
    "#2478B4",
    "#E67E22",
    "#2E9D63",
    "#B34C8C",
    "#8A6D3B",
    "#6C63B5",
    "#C44E52",
    "#4C9A9A",
]
SHEET_NAMES = [
    "Fibronectin Plot",
    "Alignment Plot",
    "Area Plot",
    "StdDev Plot",
    "Min Plot",
    "Max Plot",
    "Median Plot",
    "Alignment Filtered",
    "Area Filtered",
    "StdDev Filtered",
    "Min Filtered",
    "Max Filtered",
    "Median Filtered",
    "Merged Data",
    "Filtered Data",
    "Excluded Data",
    "Filter Summary",
    "Plate Map",
    "QC",
    "Run Log",
]


class ValidationError(Exception):
    """Carry a validation message and all affected records."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = details or []


class EventLogger(Protocol):
    """
    A report stage records events without owning a logging destination.
    """

    def event(
        self,
        level: str,
        stage: str,
        message: object,
        console: bool = True,
        timestamp: str | None = None,
    ) -> object:
        """Record a stage event at the configured destination."""
        ...


class ParsedRecord(TypedDict, total=False):
    """
    One literal input row and its validated, optional source metadata.
    """

    source_row: int
    raw: dict[str, str]
    image_id: str
    well: str
    image_number: str | None
    sequence_number: str | None
    numbers: dict[str, float]


class ReportData(TypedDict):
    """
    Complete observations and derived views, without measurement
    conversion.
    """

    columns: list[str]
    rows: list[dict[str, Any]]
    metric: str
    angle_label: str
    retained_rows: list[dict[str, Any]]
    excluded_rows: list[dict[str, Any]]
    fn_threshold: float
    group_filter_counts: list[dict[str, Any]]
    well_filter_counts: list[dict[str, Any]]
    group_order: list[str]
    group_wells: dict[str, list[str]]
    well_counts: dict[str, int]
    well_map: dict[str, str]
    plate_matrix: list[list[Any]]
    plate_id: str
    template_sheet: str
    thickness_units: dict[str, str]
    qc: list[dict[str, Any]]
    field_map: list[dict[str, str]]


class ReportInputs(TypedDict):
    """Validated paths to the three summaries and literal plate map."""

    alignment: Path
    thickness: Path
    fibronectin: Path
    template: Path
