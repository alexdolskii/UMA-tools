"""Compatibility import for the moved report implementation."""

import sys

if __name__ == "__main__":
    from .cli import report

    raise SystemExit(report())
else:
    from .reporting import workflow as _implementation

    sys.modules[__name__] = _implementation
