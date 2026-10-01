"""Core output location and file operations with explicit policies."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any


def assay_directory(source: Path, *, create: bool = False) -> Path:
    """Locate the core workflows' only results root below original images.

    Never create a missing source or follow a linked output container.
    Readers use this same location without searching historical layouts.
    """
    source = Path(source)
    if not source.is_dir():
        raise FileNotFoundError(f"Source folder not found: {source}")
    directory = source / "uma_assay"
    if directory.is_symlink():
        raise OSError(
            f"Linked uma_assay directory is not accepted: {directory}"
        )
    if create:
        directory.mkdir(exist_ok=True)
    if not directory.is_dir():
        raise FileNotFoundError(
            f"UMA results directory not found: {directory}"
        )
    return directory


def save_json(
    path: Path,
    data: Any,
    *,
    atomic: bool = True,
    allow_nan: bool = True,
    trailing_newline: bool = True,
    temporary_suffix: str = ".tmp",
    encoding: str = "utf-8",
) -> None:
    """Save JSON in each workflow's established byte format.

    Defaults reproduce the collector format. Area uses ``.pending``, no
    final newline and rejects nonfinite numbers. Report tables use the
    same strict JSON representation with ``atomic=False``.
    """
    target = path.with_name(path.name + temporary_suffix) if atomic else path
    text = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=allow_nan)
    if trailing_newline:
        text += "\n"
    target.write_text(text, encoding=encoding)
    if atomic:
        target.replace(path)


def save_csv(
    path: Path,
    columns: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
    *,
    encoding: str = "utf-8",
) -> None:
    """Write ordered records with the requested encoding."""
    with path.open("w", encoding=encoding, newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    """Hash a file in bounded chunks without changing its contents."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_label(name: str) -> str:
    """Keep readable folder names within portable filename lengths."""
    label = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", name).strip(" .")
    label = label or "images"
    encoded = label.encode("utf-8")
    if len(encoded) > 110:
        suffix = hashlib.sha256(encoded).hexdigest()[:8]
        label = encoded[:100].decode("utf-8", errors="ignore") + "_" + suffix
    return label
