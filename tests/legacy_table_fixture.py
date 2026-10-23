"""Isolate existing CSV/math regressions from the new scheme/route boundary.

These fixtures intentionally contain historical table-only inputs. The real
five-command boundary and provenance are exercised in test_shared_plate.py.
This helper never patches measurement, collection, filtering or statistics.
"""

from contextlib import ExitStack, nullcontext
from types import SimpleNamespace
from unittest.mock import patch


def table_only_boundary(test_case):
    from uma_tools.plate_templates import preflight_templates

    def preflight(folders, args, validator=None):
        if validator is None:
            prepared = [
                SimpleNamespace(source=p, error=None, selection=None)
                for p in folders
            ]
        else:
            prepared = preflight_templates(
                folders, validator, interactive=False
            )
        return [(item, None) for item in prepared]

    stack = ExitStack()
    test_case._table_boundary = stack
    test_case.addCleanup(stack.close)
    stack.enter_context(patch("uma_tools.fiber_routes.preflight", preflight))
    stack.enter_context(
        patch("uma_tools.fiber_routes.routes", return_value=[None])
    )
    stack.enter_context(
        patch(
            "uma_tools.fiber_routes.use_route",
            side_effect=lambda _: nullcontext(),
        )
    )
