#!/usr/bin/env python3
"""Compatibility launcher; prefer 2_thickness.py."""

import sys

if __name__ == "__main__":
    from uma_tools.cli import thickness

    raise SystemExit(thickness())
else:
    from uma_tools import thickness_analysis as _implementation

    sys.modules[__name__] = _implementation
