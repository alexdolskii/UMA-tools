"""Legacy import adapter for report helpers."""

import sys

from .reporting import engine as _implementation

sys.modules[__name__] = _implementation
