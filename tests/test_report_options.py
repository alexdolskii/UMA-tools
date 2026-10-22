"""Guided setup and unattended reporting have identical option semantics."""

import argparse
import sys

import pytest

from uma_tools import report
from uma_tools.report_options import parse_args, resolve_options


def terminal(monkeypatch, answers):
    prompts = []
    values = iter(answers)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)

    def answer(prompt):
        prompts.append(prompt)
        return next(values)

    monkeypatch.setattr("builtins.input", answer)
    return prompts


def test_bare_command_defaults_and_no_double_prompt(tmp_path, monkeypatch):
    path = tmp_path / "input_paths.json"
    path.write_text('{"folder_paths": ["' + str(tmp_path) + '"]}')
    prompts = terminal(monkeypatch, [str(path), "", "", "", "", "", "", ""])
    args = parse_args([])
    assert not prompts  # Supervisor preflight is deliberately read-only.
    resolved = resolve_options(args)
    assert len(prompts) == 8
    assert resolved.input == str(path)
    assert (
        resolved.fn_threshold,
        resolved.with_fia_nuclei,
        resolved.min_nuclei,
    ) == (20, False, 0)
    assert resolved.stats_unit is None
    assert resolved.plot_format == "pdf"


def test_fia_and_image_statistics_can_be_selected_by_dialog(monkeypatch):
    prompts = terminal(
        monkeypatch,
        [
            "abc",
            "101",
            "30",
            "yes",
            "-1",
            "2.5",
            "20",
            "image",
            "both",
            "Plate",
            "Experiment",
            "y",
        ],
    )
    args = resolve_options(parse_args(["-i", "inputs.json"]))
    assert args.with_fia_nuclei and args.min_nuclei == 20
    assert args.fn_threshold == 30 and args.stats_unit == "image"
    assert args.plot_format == "both" and args.sheet == "Plate"
    assert args.plate_id == "Experiment"
    assert len(prompts) == 12


def test_explicit_settings_are_never_overwritten(monkeypatch):
    prompts = terminal(monkeypatch, ["", "", "y"])
    args = resolve_options(
        parse_args(
            [
                "-i",
                "input.json",
                "--interactive",
                "--fn-threshold",
                "10",
                "--with-fia-nuclei",
                "--min-nuclei",
                "7",
                "--stats-unit",
                "well",
                "--plot-format",
                "png",
            ]
        )
    )
    assert len(prompts) == 3  # worksheet, plate label, confirmation only
    assert args.fn_threshold == 10 and args.min_nuclei == 7
    assert args.stats_unit == "well"


@pytest.mark.parametrize(
    "arguments",
    [
        ["-i", "input.json", "--non-interactive"],
        ["-i", "input.json", "--fn-threshold", "20"],
    ],
)
def test_cli_never_prompts_when_options_select_unattended_mode(
    monkeypatch, arguments
):
    terminal(monkeypatch, [])
    args = resolve_options(parse_args(arguments))
    assert not args.wizard


@pytest.mark.parametrize(
    "arguments",
    [
        [],
        ["--interactive"],
        ["-i", "x", "--min-nuclei", "20"],
        ["-i", "x", "--fn-threshold", "nan"],
        ["-i", "x", "--min-nuclei", "-1"],
    ],
)
def test_nonterminal_errors_are_actionable(monkeypatch, arguments):
    monkeypatch.setattr(sys.stdin, "isatty", lambda: False)
    with pytest.raises(SystemExit) as error:
        parse_args(arguments)
    assert error.value.code == 2


@pytest.mark.parametrize(
    "answer", ["q", "no", EOFError(), KeyboardInterrupt()]
)
def test_cancel_before_output_creation(monkeypatch, answer):
    args = argparse.Namespace(
        wizard=True,
        input="x",
        supplied={
            "--fn-threshold",
            "--with-fia-nuclei",
            "--min-nuclei",
            "--stats-unit",
            "--plot-format",
            "--sheet",
            "--plate-id",
        },
        fn_threshold=20,
        with_fia_nuclei=False,
        min_nuclei=0,
        stats_unit=None,
        plot_format="pdf",
    )
    monkeypatch.setattr(report, "parse_args", lambda argv: args)

    def reply(prompt):
        if isinstance(answer, BaseException):
            raise answer
        return answer

    monkeypatch.setattr("builtins.input", reply)
    monkeypatch.setattr(
        report,
        "read_config",
        lambda path: pytest.fail("Cancellation reached processing"),
    )
    assert report.main([]) == 130
