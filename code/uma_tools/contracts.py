"""Stable names and units used at the boundaries between UMA stages."""

IMAGE_EXTENSIONS = frozenset({".nd2", ".tif", ".tiff"})
ALIGNMENT_SUFFIX = "_processed_orientation_distribution.csv"
THICKNESS_METRICS = ("Area", "StdDev", "Min", "Max", "Median")
THICKNESS_UNITS = {
    "Area": "µm²",
    "StdDev": "µm",
    "Min": "µm",
    "Max": "µm",
    "Median": "µm",
}
EVENT_COLUMNS = ["Timestamp_UTC", "Level", "Stage", "Message"]
SUMMARY_NAMES = {
    "Alignment": "Alignment_Summary.csv",
    "Thickness": "Thickness_Summary.csv",
    "Area": "Fibronectin_Area_Summary.csv",
}

# Legacy directories ended at seconds. New image-analysis runs append
# microseconds and the process ID to avoid overwriting a prior run.
ASSAY_TIMESTAMP_PATTERN = r"(\d{8}_\d{6})(?:_(\d{6})_\d+(?:_\d+)?)?"
