"""Legacy import adapter for alignment."""

import sys

if __name__ == "__main__":
    from .cli import alignment

    raise SystemExit(alignment())
else:
    from .assays import alignment as _implementation

    sys.modules[__name__] = _implementation
