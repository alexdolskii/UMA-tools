#!/usr/bin/env python3
"""Compatibility launcher; prefer 4_collect_results.py."""

import sys

if __name__ == "__main__":
    from uma_tools.cli import collect_results

    raise SystemExit(collect_results())
else:
    from uma_tools import collect_results as _implementation

    sys.modules[__name__] = _implementation
