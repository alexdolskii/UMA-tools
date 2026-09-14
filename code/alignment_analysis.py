#!/usr/bin/env python3
"""Compatibility launcher; prefer 1_alignment.py."""

import sys

if __name__ == "__main__":
    from uma_tools.cli import alignment

    raise SystemExit(alignment())
else:
    from uma_tools import alignment_analysis as _implementation

    sys.modules[__name__] = _implementation
