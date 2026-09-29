#!/usr/bin/env python3
"""Compatibility launcher for the installed uma_nla_segment command."""

from pathlib import Path

from nuclei_layers_assay.segment import main

if __name__ == "__main__":
    raise SystemExit(
        main(default_input=str(Path(__file__).with_name("nuclei_layers.json")))
    )
