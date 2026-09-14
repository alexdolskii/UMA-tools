"""Legacy import adapter for thickness."""

import sys

if __name__ == "__main__":
    from .cli import thickness

    raise SystemExit(thickness())
else:
    from .assays import thickness as _implementation

    sys.modules[__name__] = _implementation
