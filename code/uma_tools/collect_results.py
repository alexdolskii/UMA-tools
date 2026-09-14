"""Compatibility import for the moved collect_results implementation."""

import sys

if __name__ == "__main__":
    from .cli import collect_results

    raise SystemExit(collect_results())
else:
    from . import collection as _implementation

    sys.modules[__name__] = _implementation
