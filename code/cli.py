"""Compatibility import for UMA command entry points."""

import sys

from uma_tools import cli as _implementation

sys.modules[__name__] = _implementation
