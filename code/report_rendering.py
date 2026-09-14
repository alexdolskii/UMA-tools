"""Compatibility import for reusable report helpers."""

import sys

from uma_tools import report_rendering as _implementation

sys.modules[__name__] = _implementation
