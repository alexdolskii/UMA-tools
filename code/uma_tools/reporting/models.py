"""Public data contracts for report validation and rendering."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, TypedDict


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
