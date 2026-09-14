#!/usr/bin/env python3
"""Compatibility launcher; prefer 5_report.py."""

import sys

if __name__ == "__main__":
    from uma_tools.cli import report

    raise SystemExit(report())
else:
    from uma_tools import report as _implementation

    sys.modules[__name__] = _implementation
