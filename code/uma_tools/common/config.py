"""Configuration primitives that preserve each assay's path policies."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from .errors import ValidationError


def reject_metadata_json(
    path: Path,
    *,
    error_type: type[Exception] = ValidationError,
    message: str | None = None,
) -> None:
    """Reject AppleDouble configurations before opening them."""
    if path.name.startswith("._"):
        raise error_type(
            message
            or "macOS metadata JSON files (._...) are not input configurations"
        )


def load_json(
    path: Path,
    *,
    encoding: str = "utf-8-sig",
    reject_metadata: bool = True,
    error_type: type[Exception] = ValidationError,
) -> Any:
    """Read JSON with explicit encoding and metadata rejection."""
    if reject_metadata:
        reject_metadata_json(path, error_type=error_type)
    return json.loads(path.read_text(encoding=encoding))


def folder_entries(
    data: Any,
    *,
    error_type: type[Exception] = ValidationError,
    list_message: str = "JSON must contain a nonempty folder_paths list",
    entry_message: str = (
        "Each folder_paths entry must be a nonempty path string"
    ),
) -> list[str]:
    """Validate entries, preserving relative paths and duplicates."""
    folders = data.get("folder_paths") if isinstance(data, dict) else None
    if not isinstance(folders, list) or not folders:
        raise error_type(list_message)
    if any(
        not isinstance(folder, str) or not folder.strip() for folder in folders
    ):
        raise error_type(entry_message)
    return folders


def resolve_path(value: str, relative_to: Path) -> Path:
    """Resolve a path relative to an explicit working directory."""
    path = Path(value).expanduser()
    return (path if path.is_absolute() else relative_to / path).resolve()


def resolve_folders(
    values: list[str],
    *,
    relative_to: Path | None = None,
    use_realpath: bool = False,
) -> list[Path]:
    """Resolve source entries, preserving order and repeats.

    ``relative_to=None`` follows the current working directory. The
    collector does not resolve symbolic links until duplicate detection;
    Area explicitly requests ``use_realpath=True`` relative to its JSON.
    """
    folders = []
    for value in values:
        expanded = os.path.expanduser(value)
        if relative_to is not None and not os.path.isabs(expanded):
            expanded = os.path.join(relative_to, expanded)
        path = Path(os.path.abspath(expanded))
        folders.append(path.resolve() if use_realpath else path)
    return folders


def read_config(path: Path) -> list[Path]:
    """Read collector/report configuration using CWD-relative paths."""
    reject_metadata_json(path)
    if path.suffix.lower() != ".json":
        raise ValidationError("The input configuration must be a JSON file")
    data = load_json(path, reject_metadata=False)
    return resolve_folders(folder_entries(data))
