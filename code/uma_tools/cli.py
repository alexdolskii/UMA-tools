"""Console commands with lazy scientific imports."""

import argparse
import logging
import sys
from collections.abc import Callable, Sequence
from typing import Any

from .common.version import package_version


def _shutdown_imagej_workers() -> None:
    """Delegate legacy worker cleanup to runtime infrastructure."""
    from .runtime.imagej import shutdown_imagej_workers

    shutdown_imagej_workers()


def _run_imagej_command(callback: Callable[..., Any], *args: Any) -> int:
    """Finish workers without masking the original analysis error."""
    try:
        callback(*args)
        return 0
    finally:
        analysis_failed = sys.exc_info()[0] is not None
        try:
            _shutdown_imagej_workers()
        except Exception:
            if not analysis_failed:
                raise
            logging.exception(
                "Could not close ImageJ workers after analysis failed."
            )


def _image_parser(description: str) -> argparse.ArgumentParser:
    """Build shared input/version options without starting Fiji."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {package_version()}",
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to a JSON file containing folder_paths",
    )
    return parser


def alignment(argv: Sequence[str] | None = None) -> int:
    """Run alignment with its existing angle and prompts."""
    parser = _image_parser("Fibronectin alignment analysis")
    parser.add_argument(
        "-a",
        "--angle_value",
        type=float,
        default=15,
        help="Alignment angle in degrees (default: 15)",
    )
    args = parser.parse_args(argv)
    from .assays.alignment import main_fibronectin_processing

    return _run_imagej_command(
        main_fibronectin_processing, args.input, args.angle_value
    )


def thickness(argv: Sequence[str] | None = None) -> int:
    """Run thickness and finish its ImageJ workers."""
    parser = _image_parser("Fibronectin thickness analysis")
    args = parser.parse_args(argv)
    print(f"UMA-tools {package_version()} — thickness analysis", flush=True)
    from .assays.thickness import main

    return _run_imagej_command(main, args.input)


def area(argv: Sequence[str] | None = None) -> int:
    """Run area and record its runtime shutdown status."""
    from .assays.area import main

    return main() if argv is None else main(argv)


def collect_results(argv: Sequence[str] | None = None) -> int:
    """Collect summaries without loading scientific runtimes."""
    from .collection import main

    return main() if argv is None else main(argv)


def report(argv: Sequence[str] | None = None) -> int:
    """Create plots and an Excel report without starting Fiji."""
    from .reporting.workflow import main

    return main() if argv is None else main(argv)
