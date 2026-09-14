#!/usr/bin/env python3
"""Compatibility launcher; prefer 3_area.py."""

import sys

if __name__ == "__main__":
    from uma_tools.cli import area

    raise SystemExit(area())
else:
    from uma_tools import area_analysis as _implementation

    sys.modules[__name__] = _implementation
