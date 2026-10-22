"""Report options and a terminal wizard without scientific imports."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from . import package_version
from .config import read_config


def nonnegative_integer(value):
    """Accept a literal nonnegative integer, including zero."""
    try:
        result = int(value)
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError(
            "Enter a nonnegative integer"
        ) from error
    if result < 0:
        raise argparse.ArgumentTypeError("Enter a nonnegative integer")
    return result


def coverage(value):
    result = float(value)
    if not math.isfinite(result) or not 0 <= result <= 100:
        raise argparse.ArgumentTypeError(
            "FN coverage must be between 0 and 100"
        )
    return result


def parse_args(argv=None):
    """Parse only: the runtime supervisor must never prompt twice."""
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Create UMA plots and Excel reports; "
            "run without options for guided setup"
        ),
    )
    parser.add_argument("-i", "--input", help="JSON containing folder_paths")
    parser.add_argument(
        "--fn-threshold",
        type=coverage,
        default=20.0,
        help="Minimum FN coverage, percent (default: 20; inclusive)",
    )
    parser.add_argument(
        "--with-fia-nuclei",
        action="store_true",
        help="Use the latest valid final-nuclei run in each source/fia_assay",
    )
    parser.add_argument(
        "--min-nuclei",
        type=nonnegative_integer,
        default=0,
        help=(
            "Minimum non-border nuclei per image "
            "(default: 0/off; requires FIA)"
        ),
    )
    parser.add_argument(
        "--stats-unit",
        choices=("well", "image"),
        default=None,
        help="Optional Welch/Holm tests after filtering; omitted = no tests",
    )
    parser.add_argument(
        "--plot-format", choices=("pdf", "png", "both"), default="pdf"
    )
    parser.add_argument(
        "--sheet", default=None, help="Plate worksheet (default: first)"
    )
    parser.add_argument(
        "--plate-id", default="", help="Label (default: source folder name)"
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--interactive",
        action="store_true",
        help="Ask for unspecified settings",
    )
    mode.add_argument(
        "--non-interactive",
        action="store_true",
        help="Use defaults without questions",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {package_version()}"
    )
    tokens = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(tokens)
    args.supplied = {
        (
            "-i"
            if token.startswith("-i") and not token.startswith("--")
            else token.split("=", 1)[0]
        )
        for token in tokens
        if token.startswith("-")
    }
    only_input = args.supplied <= {"-i", "--input"}
    args.wizard = args.interactive or (
        sys.stdin.isatty() and not args.non_interactive and only_input
    )
    if args.wizard and not sys.stdin.isatty():
        parser.error(
            "Interactive setup needs a terminal; "
            "supply options and --non-interactive"
        )
    if not args.wizard and not args.input:
        parser.error("-i/--input is required outside interactive setup")
    if args.min_nuclei and not args.with_fia_nuclei and not args.wizard:
        parser.error("--min-nuclei requires --with-fia-nuclei")
    return args


def _ask(label, default, convert=str):
    """Retry invalid answers; q, EOF and Ctrl-C cancel before processing."""
    while True:
        answer = input(f"{label} [{default}]: ").strip()
        if answer.lower() in ("q", "quit"):
            raise KeyboardInterrupt
        try:
            return convert(answer if answer else str(default))
        except (ValueError, OSError, argparse.ArgumentTypeError) as error:
            print(f"Invalid answer: {error}. Please try again (q = cancel).")


def _choice(options):
    def convert(value):
        key = value.lower()
        if key not in options:
            raise ValueError("Choose " + ", ".join(options))
        return options[key]

    return convert


def _config_path(value):
    path = Path(value.strip("\"'")).expanduser()
    read_config(path)
    return str(path)


def resolve_options(args):
    """Ask once in the worker and return the resolved command options."""
    if not args.wizard:
        return args
    print("UMA report setup. Enter accepts the default; q cancels.")
    if not args.input:
        args.input = _ask("Input JSON path", "input_paths.json", _config_path)
    if "--fn-threshold" not in args.supplied:
        args.fn_threshold = _ask("Minimum FN coverage (%)", 20, coverage)
    if "--with-fia-nuclei" not in args.supplied:
        args.with_fia_nuclei = _ask(
            "Use FIA nuclei counts from the same source folders? y/n",
            "y" if args.min_nuclei else "n",
            _choice({"y": True, "yes": True, "n": False, "no": False}),
        )
    if args.min_nuclei and not args.with_fia_nuclei:
        raise ValueError("--min-nuclei requires FIA nuclei integration")
    if args.with_fia_nuclei and "--min-nuclei" not in args.supplied:
        args.min_nuclei = _ask(
            "Minimum non-border nuclei per image (0 = off)",
            0,
            nonnegative_integer,
        )
    if "--stats-unit" not in args.supplied:
        args.stats_unit = _ask(
            "Statistics: none, well, or image",
            "none",
            _choice({"none": None, "well": "well", "image": "image"}),
        )
    if "--plot-format" not in args.supplied:
        args.plot_format = _ask(
            "Figure files: pdf, png, or both",
            "pdf",
            _choice({key: key for key in ("pdf", "png", "both")}),
        )
    print(
        f"Input: {args.input}\nFN coverage >= {args.fn_threshold:g}%; "
        f"FIA nuclei: {'enabled' if args.with_fia_nuclei else 'disabled'}; "
        f"minimum nuclei: {args.min_nuclei}; "
        f"statistics: {args.stats_unit or 'disabled'}; "
        f"figures: {args.plot_format}"
    )
    if not _ask(
        "Create reports? y/n",
        "y",
        _choice({"y": True, "yes": True, "n": False, "no": False}),
    ):
        raise KeyboardInterrupt
    return args
