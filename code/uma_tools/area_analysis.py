"""Compatibility import for the moved area_analysis implementation."""

import sys

if __name__ == "__main__":
    from .cli import area

    raise SystemExit(area())
else:
    from .assays import area as _implementation

    sys.modules[__name__] = _implementation
