"""Select and validate source plate maps before a report batch starts."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .files import sha256_file
from .progress import console, prompt_value
from .report_schema import ValidationError

Validator = Callable[[Path, Path], str]


def candidates(source: Path, prefix="UMA") -> list[Path]:
    """List only direct, regular UMA*.xlsx inputs; never search results."""
    if not source.is_dir():
        raise ValidationError(f"Source folder not found: {source}")
    return sorted(
        (
            path
            for path in source.iterdir()
            if path.name.casefold().startswith(prefix.casefold())
            and path.suffix.casefold() == ".xlsx"
            and not path.is_symlink()
            and path.is_file()
        ),
        key=lambda path: (path.name.casefold(), path.name),
    )


def find_template(source: Path, prefix="UMA") -> Path:
    """Resolve an unambiguous map without asking or inspecting nested runs."""
    paths = candidates(source, prefix)
    if len(paths) != 1:
        listing = "; ".join(str(path) for path in paths) or "none"
        raise ValidationError(
            f"Expected one {prefix}*.xlsx plate template directly in "
            f"{source}; "
            f"found {len(paths)}: {listing}. "
            "Use an interactive terminal to select a template, or leave "
            "one matching file in this source folder."
        )
    return paths[0]


@dataclass(frozen=True)
class TemplateSelection:
    source: Path
    path: Path
    sha256: str
    worksheet: str
    method: str

    def verify(self) -> None:
        """Do not let different variants use different edits of a plate map."""
        if (
            self.path.is_symlink()
            or not self.path.is_file()
            or sha256_file(self.path) != self.sha256
        ):
            raise ValidationError(
                f"Plate template changed after preflight: {self.path}. "
                "Rerun the report to validate the edited template."
            )

    def verify_snapshot(self, snapshot: Path) -> None:
        if sha256_file(snapshot) != self.sha256:
            raise ValidationError(
                f"Archived plate template differs from preflight: {self.path}"
            )
        self.verify()

    def metadata(self) -> dict:
        return {
            "source_folder": str(self.source),
            "original_path": str(self.path),
            "sha256": self.sha256,
            "worksheet": self.worksheet,
            "selection_method": self.method,
        }


@dataclass(frozen=True)
class PreflightResult:
    source: Path
    selection: TemplateSelection | None = None
    error: ValidationError | None = None


def prepare_template(
    source: Path,
    validator: Validator,
    path: Path | None = None,
    method: str = "automatic",
    prefix: str = "UMA",
) -> TemplateSelection:
    """Validate content and capture the exact bytes approved at preflight."""
    path = find_template(source, prefix) if path is None else path
    if path.parent != source or path not in candidates(source, prefix):
        raise ValidationError(
            f"Not a regular source-folder {prefix}*.xlsx: {path}"
        )
    digest = sha256_file(path)
    # Workbook parsing is an input boundary: corrupt ZIP/XML, unavailable
    # worksheets and design errors must be fixable before any report starts.
    try:
        sheet = validator(path, source)
    except Exception as error:
        details = getattr(error, "details", None)
        description = str(error)
        if details:
            description += " " + "; ".join(
                json.dumps(detail, ensure_ascii=False) for detail in details
            )
        raise ValidationError(
            f"Invalid plate template {path}: {description}", details
        ) from error
    selection = TemplateSelection(source, path, digest, sheet, method)
    selection.verify()
    return selection


def interactive_selection(args) -> bool:
    """Template questions are independent of the analysis-settings wizard."""
    return sys.stdin.isatty() and not getattr(args, "non_interactive", False)


def _answer(message: str) -> str:
    value = prompt_value(message).strip().casefold()
    if value in ("q", "quit", "cancel"):
        raise KeyboardInterrupt
    return value


def preflight_templates(
    folders: list[Path],
    validator: Validator,
    *,
    interactive: bool,
    prefix="UMA",
) -> list[PreflightResult]:
    """Resolve every source before returning control to report generation."""
    results, seen = [], set()
    console("Plate-template preflight: checking all source folders.")
    for source in folders:
        canonical = source.resolve()
        if canonical in seen:
            console(f"Skipping duplicate JSON folder: {source}")
            continue
        seen.add(canonical)
        while True:
            try:
                paths = candidates(source, prefix)
                path, method = None, "automatic"
                if len(paths) > 1 and interactive:
                    console(f"Several {prefix} plate templates in {source}:")
                    for index, candidate in enumerate(paths, 1):
                        console(f"  {index}. {candidate.name}")
                    answer = _answer(
                        "Template number, r = rescan, s = skip folder, "
                        "q = cancel: "
                    )
                    if answer == "r":
                        continue
                    if answer == "s":
                        issue = ValidationError(
                            f"Folder skipped: ambiguous plate templates in "
                            f"{source}: " + ", ".join(p.name for p in paths)
                        )
                        results.append(PreflightResult(source, error=issue))
                        console(str(issue))
                        break
                    if not answer.isdecimal() or not 1 <= int(answer) <= len(
                        paths
                    ):
                        console("Enter a template number from the list.")
                        continue
                    path, method = paths[int(answer) - 1], "interactive"
                selection = prepare_template(
                    source, validator, path, method, prefix
                )
                results.append(PreflightResult(source, selection))
                console(f"READY: {source} -> {selection.path.name}")
                break
            except (OSError, ValidationError) as error:
                issue = (
                    error
                    if isinstance(error, ValidationError)
                    else ValidationError(str(error))
                )
                console(f"NOT READY: {source}: {issue}")
                if interactive:
                    answer = _answer(
                        "Fix this folder, then r = retry [default], "
                        "s = skip folder, q = cancel: "
                    )
                    if answer in ("", "r", "retry"):
                        continue
                    if answer not in ("s", "skip"):
                        console("Enter r, s or q.")
                        continue
                results.append(PreflightResult(source, error=issue))
                break
    console("Plate-template preflight summary (before processing):")
    for result in results:
        detail = (
            f"READY -> {result.selection.path}"
            if result.selection
            else f"SKIPPED / NOT READY -> {result.error}"
        )
        console(f"  {result.source}: {detail}")
    return results
