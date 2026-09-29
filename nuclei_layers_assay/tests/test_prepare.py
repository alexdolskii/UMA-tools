"""Behavioral tests without Java or the later-stage dependencies."""

import json
import logging
import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest
import uma_tools.cli
from nuclei_layers_assay import prepare


def configuration(tmp_path, folders, **options):
    path = tmp_path / "settings.json"
    path.write_text(
        json.dumps({"folder_paths": list(map(str, folders)), **options})
    )
    return path


def source_folder(tmp_path, name):
    folder = tmp_path / name
    folder.mkdir()
    (folder / "image.tif").write_bytes(b"test fixture")
    return folder


def test_defaults_and_shared_json_keep_later_settings_optional(tmp_path):
    folder = source_folder(tmp_path, "source")
    path = configuration(tmp_path, [folder, folder], model_path="unused")
    settings = prepare.read_settings(path)
    assert settings == prepare.Settings((folder,), 1, 4.0, 3)


@pytest.mark.parametrize(
    "options",
    [
        {"nuclei_channel": 0},
        {"nuclei_channel": True},
        {"nuclei_channel": 1.5},
        {"gaussian_sigma": float("nan")},
        {"gaussian_sigma": -1},
        {"mean_radius": -1},
        {"mean_radius": 2.5},
    ],
)
def test_invalid_filter_or_channel_is_rejected(tmp_path, options):
    path = configuration(tmp_path, [tmp_path], **options)
    with pytest.raises(ValueError):
        prepare.read_settings(path)


def test_metadata_json_rejected_before_opening(tmp_path):
    with pytest.raises(ValueError, match="metadata JSON"):
        prepare.read_settings(tmp_path / "._missing.json")


def test_missing_folder_cannot_silently_drop_an_experiment(tmp_path):
    path = configuration(tmp_path, [tmp_path, tmp_path / "missing"])
    with pytest.raises(ValueError, match="does not exist"):
        prepare.read_settings(path)


def test_discovery_excludes_metadata_hidden_files_and_directories(tmp_path):
    folder = source_folder(tmp_path, "source")
    for name in ("._image.tif", ".hidden.nd2", "_other.tif", "notes.txt"):
        (folder / name).touch()
    (folder / "directory.tif").mkdir()
    (folder / "second.ND2").touch()
    settings = prepare.Settings((folder,), 1, 4.0, 3)
    assert [p.name for p in prepare.discover_inputs(settings)[folder]] == [
        "image.tif",
        "second.ND2",
    ]


def test_output_name_collision_rejected_before_overwrite(tmp_path):
    folder = source_folder(tmp_path, "source")
    (folder / "IMAGE.nd2").touch()
    settings = prepare.Settings((folder,), 1, 4.0, 3)
    with pytest.raises(ValueError, match="same nuclei TIFF"):
        prepare.discover_inputs(settings)


def test_source_inside_output_is_never_deleted(tmp_path):
    folder = source_folder(tmp_path, "source")
    nested = source_folder(folder, "processed")
    settings = prepare.Settings((folder, nested), 1, 4.0, 3)
    with pytest.raises(ValueError, match="contains a source"):
        prepare.discover_inputs(settings)
    assert (nested / "image.tif").exists()


def test_no_replacement_without_confirmation(tmp_path, monkeypatch):
    folder = source_folder(tmp_path, "source")
    (folder / "processed").mkdir()
    old = folder / "processed" / "old.tif"
    old.write_bytes(b"previous result")
    path = configuration(tmp_path, [folder])
    initialize = Mock(side_effect=AssertionError("must not start Fiji"))
    monkeypatch.setattr(prepare, "initialize_imagej", initialize)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: False)
    assert prepare.main(["-i", str(path)]) == 1
    assert old.read_bytes() == b"previous result"
    initialize.assert_not_called()


def test_one_prompt_for_all_folders(tmp_path, monkeypatch):
    folders = [source_folder(tmp_path, name) for name in ("first", "second")]
    for folder in folders:
        (folder / "processed").mkdir()
    prompt = Mock(return_value="yes")
    monkeypatch.setattr("builtins.input", prompt)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: True)
    prepare.confirm_replacement(
        prepare.Settings(tuple(folders), 1, 4, 3), overwrite=False
    )
    prompt.assert_called_once()


def test_fiji_startup_failure_preserves_old_results(tmp_path, monkeypatch):
    folder = source_folder(tmp_path, "source")
    (folder / "processed").mkdir()
    old = folder / "processed" / "old.tif"
    old.write_bytes(b"previous result")
    path = configuration(tmp_path, [folder])
    monkeypatch.setattr(
        prepare,
        "initialize_imagej",
        Mock(side_effect=RuntimeError("JVM unavailable")),
    )
    shutdown = Mock()
    monkeypatch.setattr(uma_tools.cli, "_shutdown_imagej_workers", shutdown)
    assert prepare.main(["-i", str(path), "--overwrite"]) == 1
    assert old.read_bytes() == b"previous result"
    shutdown.assert_called_once()


@pytest.mark.parametrize("image_failure", [False, True])
def test_all_folders_finish_and_workers_stop_once(
    tmp_path,
    monkeypatch,
    image_failure,
):
    folders = [source_folder(tmp_path, name) for name in ("first", "second")]
    path = configuration(tmp_path, folders)
    context = Mock()
    initialize = Mock(return_value=context)
    monkeypatch.setattr(prepare, "initialize_imagej", initialize)
    visited = []

    class Processor:
        def process(self, source, target, settings):
            visited.append(source.parent)
            target.write_bytes(b"output")
            if image_failure and source.parent == folders[0]:
                raise ValueError("unreadable image")
            return "one complete stack"

    monkeypatch.setattr(prepare, "NucleiProcessor", Processor)
    shutdown = Mock()
    monkeypatch.setattr(uma_tools.cli, "_shutdown_imagej_workers", shutdown)
    handlers = list(logging.getLogger("nuclei_layers_assay.prepare").handlers)
    assert prepare.main(["-i", str(path)]) == int(image_failure)
    assert visited == folders
    initialize.assert_called_once()
    context.dispose.assert_called_once()
    shutdown.assert_called_once()
    assert (folders[1] / "processed" / "image_nuclei.tif").exists()
    if image_failure:
        assert not (folders[0] / "processed" / "image_nuclei.tif").exists()
    log = (folders[0] / "nuclei_analysis1.log").read_text()
    assert str(folders[1]) not in log
    assert (
        logging.getLogger("nuclei_layers_assay.prepare").handlers == handlers
    )


@pytest.mark.parametrize("option", ["--help", "--version"])
def test_installed_cli_outside_repository_needs_no_jvm(tmp_path, option):
    script = (
        "import sys; from nuclei_layers_assay.prepare import main; "
        f"\ntry: main([{option!r}])\nexcept SystemExit as e: "
        "assert e.code == 0\n"
        "assert 'imagej' not in sys.modules\n"
        "assert 'tensorflow' not in sys.modules\n"
        "assert 'stardist' not in sys.modules\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_legacy_launcher_remains_available(tmp_path):
    script = Path(__file__).parents[1] / "1_nla_fiji_channel_extraction.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=20,
    )
    assert result.returncode == 0, result.stdout + result.stderr
