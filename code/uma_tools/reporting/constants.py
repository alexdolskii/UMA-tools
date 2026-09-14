"""
Stable measurement names, units, identity patterns, and output schema.
"""

import re

from ..common.contracts import ALIGNMENT_SUFFIX
from ..common.contracts import EVENT_COLUMNS as EVENT_COLUMNS
from ..common.contracts import THICKNESS_METRICS as THICKNESS_METRICS
from ..common.contracts import THICKNESS_UNITS as THICKNESS_UNITS

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
