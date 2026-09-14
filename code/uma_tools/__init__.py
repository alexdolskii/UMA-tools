"""UMA image assays, validated result collection, and plate reports.

Command entry points load scientific dependencies only when needed.
"""

from importlib.metadata import PackageNotFoundError, version


def package_version() -> str:
    """Return the installed UMA version or a source-only fallback."""
    try:
        return version("uma-tools")
    except PackageNotFoundError:
        return "not installed"
