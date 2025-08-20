#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration Loader Module
Loads and validates configuration from JSON files
"""

import json
import os


def load_config(path: str) -> dict:
    """Load configuration from JSON file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    with open(path, encoding="utf-8") as fh:
        raw = json.load(fh)

    folders = [p for p in raw.get("folder_paths", []) if os.path.isdir(p)]
    if not folders:
        raise ValueError("'folder_paths' must list existing directories")

    return {
        "folders": folders,
        "nuclei_channel": int(raw.get("nuclei_channel", 1)),
        "gaussian_sigma": float(raw.get("gaussian_sigma", 4.0)),
        "mean_radius": int(raw.get("mean_radius", 3)),
        "z_scale_factor": int(raw.get("z_scale_factor", 32)),
        "n_tiles": raw.get("n_tiles", [1, 1, 1]),
        "model_path": raw.get("model_path"),
    }
