"""File operations with explicit serialization policies."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any


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
