"""Explicit dependency preflight after the report log is open."""

from __future__ import annotations

import importlib
import importlib.metadata
import json

from .models import EventLogger


def load_dependencies(log: EventLogger) -> dict[str, str]:
    """
    Load libraries and record versions without injecting module globals.

    Individual report functions also import the libraries they use
    locally. The workflow calls this preflight after creating persistent
    diagnostics.
    """
    try:
        importlib.import_module("numpy")
        matplotlib = importlib.import_module("matplotlib")
        matplotlib.use("Agg")
        importlib.import_module("matplotlib.pyplot")
        importlib.import_module("openpyxl")
        # Openpyxl uses Pillow when embedding the generated PNG images.
        importlib.import_module("PIL.Image")
    except ImportError as error:
        raise RuntimeError(
            "A report dependency is unavailable. Install the dependencies "
            "declared by this UMA release. Original error: " + str(error)
        ) from error
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "matplotlib", "openpyxl", "Pillow")
    }
    log.event("INFO", "Dependencies", json.dumps(versions))
    return versions
