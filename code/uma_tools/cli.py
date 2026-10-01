"""Console commands with lazy scientific imports."""

import argparse
import logging
import sys
from collections.abc import Callable, Sequence
from typing import Any

from . import package_version
from .progress import CANCELLATIONS, CommandSession
from .runtime import block_cleanup, managed, record_completion


def _shutdown_imagej_workers() -> None:
    """Close the shared ImageJ worker pool after a command."""
    from .imagej import shutdown_imagej_workers

    shutdown_imagej_workers()


def _run_imagej_command(callback: Callable[..., Any], *args: Any) -> int:
    """Finish workers without masking the original analysis error."""
    try:
        result = callback(*args)
        return result if isinstance(result, int) else 0
    finally:
        analysis_failed = sys.exc_info()[0] is not None
        try:
            _shutdown_imagej_workers()
        except Exception:
            block_cleanup("ImageJ worker shutdown failed")
            if not analysis_failed:
                raise
            logging.exception(
                "Could not close ImageJ workers after analysis failed."
            )


def _invoke(step, callback, *args):
    """Keep diagnostics active until ImageJ workers have finished."""
    session = CommandSession(step)
    normal = False
    try:
        with session:
            result = callback(*args)
            session.exit_code = result if isinstance(result, int) else 0
        code = session.exit_code
        normal = not session.unavailable and (
            code in (0, 130)
            or (
                code == 1
                and bool(session.folders)
                and all(
                    session.outcomes.get(str(folder)) in ("SUCCESS", "PARTIAL")
                    for folder in session.folders
                )
            )
        )
    except CANCELLATIONS:
        print("Cancelled by user.", file=sys.stderr)
        code, normal = 130, not session.unavailable
    except Exception as error:
        print(
            f"ERROR: {error}. See uma_assay/UMA_Logs for details.",
            file=sys.stderr,
        )
        code = 1
    record_completion(code, normal, session.outcomes)
    return code


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


def parse_arguments(step, argv=None):
    """Parse arguments before allocating managed runtime resources."""
    if step in ("alignment", "thickness"):
        parser = _image_parser(f"Fibronectin {step} analysis")
        if step == "alignment":
            parser.add_argument(
                "-a",
                "--angle_value",
                type=float,
                default=15,
                help="Alignment angle in degrees (default: 15)",
            )
        return parser.parse_args(argv)
    if step == "area":
        from .area_analysis import parse_args
    elif step == "collect_results":
        from .collect_results import parse_args
    elif step == "report":
        from .report import parse_args
    else:
        raise ValueError(f"Unknown command: {step}")
    return parse_args(argv)


@managed("alignment")
def alignment(argv: Sequence[str] | None = None) -> int:
    """Run alignment with its existing angle and prompts."""
    args = parse_arguments("alignment", argv)
    from .alignment_analysis import main_fibronectin_processing

    return _invoke(
        "alignment",
        _run_imagej_command,
        main_fibronectin_processing,
        args.input,
        args.angle_value,
    )


@managed("thickness")
def thickness(argv: Sequence[str] | None = None) -> int:
    """Run thickness and finish its ImageJ workers."""
    args = parse_arguments("thickness", argv)
    print(f"UMA-tools {package_version()} — thickness analysis", flush=True)
    from .thickness_analysis import main

    return _invoke("thickness", _run_imagej_command, main, args.input)


@managed("area")
def area(argv: Sequence[str] | None = None) -> int:
    """Run area and record its runtime shutdown status."""
    from .area_analysis import main

    return _invoke("area", main, *(() if argv is None else (argv,)))


@managed("collect_results")
def collect_results(argv: Sequence[str] | None = None) -> int:
    """Collect summaries without loading scientific runtimes."""
    from .collect_results import main

    return _invoke("collect_results", main, *(() if argv is None else (argv,)))


@managed("report")
def report(argv: Sequence[str] | None = None) -> int:
    """Create plots and an Excel report without starting Fiji."""
    from .report import main

    return _invoke("report", main, *(() if argv is None else (argv,)))


def diagnostics(argv=None):
    """Inspect known resources without loading scientific libraries."""
    from .diagnostics import main

    return main(argv)
