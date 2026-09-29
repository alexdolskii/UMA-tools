"""Segmentation safety, model selection and output-contract regressions."""

import json
import logging
import subprocess
import sys
from unittest.mock import Mock

import numpy as np
import pytest
import tifffile
from nuclei_layers_assay import segment, segment_models


def model_fixture(tmp_path, monkeypatch, nested=False):
    root = tmp_path / "models"
    model = root / "trained" / ("model" if nested else "")
    model.mkdir(parents=True)
    (model / "config.json").write_text(
        json.dumps({"n_dim": 3, "n_channel_in": 1, "axes": "ZYXC"})
    )
    (model / "thresholds.json").write_text('{"prob": 0.6, "nms": 0.3}')
    (model / "weights_best.h5").write_bytes(b"fake weights for unit tests")
    monkeypatch.setattr(segment_models, "MODEL_ROOT", root)
    return model


def folder_fixture(tmp_path, name="source"):
    folder = tmp_path / name
    (folder / "processed").mkdir(parents=True)
    data = np.arange(4 * 8 * 10, dtype=np.uint16).reshape(4, 8, 10)
    tifffile.imwrite(
        folder / "processed" / "sample_nuclei.tif",
        data,
        imagej=True,
        resolution=(2.5, 2),
        metadata={"axes": "ZYX", "spacing": 0.75, "unit": "um"},
    )
    return folder


def configuration(tmp_path, folders, **extra):
    path = tmp_path / "input.json"
    path.write_text(
        json.dumps({"folder_paths": list(map(str, folders)), **extra})
    )
    return path


@pytest.mark.parametrize("nested", [False, True])
def test_moved_models_resolve_by_name_or_repository_path(
    tmp_path, monkeypatch, nested
):
    model = model_fixture(tmp_path, monkeypatch, nested)
    monkeypatch.chdir(tmp_path)
    config = tmp_path / "input.json"
    for value in (
        "trained",
        "alternative_assays/stardist_models_nuclei_layers/trained",
        "nuclei_layers_assay/stardist_models_nuclei_layers/trained",
        str(model),
    ):
        assert segment_models.resolve_model(value, config) == model
    assert segment_models.validate_model(model) == {"prob": 0.6, "nms": 0.3}


def test_existing_explicit_model_wins_over_packaged_copy(
    tmp_path, monkeypatch
):
    model_fixture(tmp_path, monkeypatch)
    explicit = tmp_path / "stardist_models_nuclei_layers" / "trained"
    explicit.mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    assert (
        segment_models.resolve_model(
            "stardist_models_nuclei_layers/trained", tmp_path / "input.json"
        )
        == explicit
    )


def test_json_relative_path_and_ambiguous_relative_path(tmp_path, monkeypatch):
    model = model_fixture(tmp_path, monkeypatch)
    working = tmp_path / "working"
    working.mkdir()
    monkeypatch.chdir(working)
    config = tmp_path / "input.json"
    assert segment_models.resolve_model("models/trained", config) == model
    (working / "models" / "trained").mkdir(parents=True)
    with pytest.raises(ValueError, match="Ambiguous"):
        segment_models.resolve_model("models/trained", config)


def test_missing_weights_and_multichannel_models_are_rejected(
    tmp_path, monkeypatch
):
    model = model_fixture(tmp_path, monkeypatch)
    (model / "weights_best.h5").unlink()
    with pytest.raises(ValueError, match="weights_best"):
        segment_models.validate_model(model)
    (model / "weights_best.h5").write_bytes(b"weights")
    (model / "config.json").write_text(
        '{"n_dim": 3, "n_channel_in": 2, "axes": "ZYXC"}'
    )
    with pytest.raises(ValueError, match="single-channel"):
        segment_models.validate_model(model)


def test_unspecified_model_never_downloads_a_demo(tmp_path, monkeypatch):
    model = model_fixture(tmp_path, monkeypatch)
    config = tmp_path / "input.json"
    with pytest.raises(ValueError, match="Specify model_path"):
        segment_models.choose_model(None, config, interactive=False)
    prompt = Mock(return_value="1")
    monkeypatch.setattr("builtins.input", prompt)
    selected, _ = segment_models.choose_model(None, config, interactive=True)
    assert selected == model
    prompt.assert_called_once()


def test_discovery_filters_metadata_and_rejects_colliding_names(tmp_path):
    folder = folder_fixture(tmp_path)
    processed = folder / "processed"
    for name in ("._sample.tif", "_hidden.tif", "notes.txt"):
        (processed / name).touch()
    (processed / "directory.tif").mkdir()
    config = configuration(tmp_path, [folder, folder])
    settings, inputs = segment.read_settings(config)
    assert settings["n_tiles"] == (1, 1, 1)
    assert [p.name for p in inputs[folder]] == ["sample_nuclei.tif"]
    (processed / "SAMPLE_NUCLEI.tiff").touch()
    with pytest.raises(ValueError, match="same mask"):
        segment.read_settings(config)


@pytest.mark.parametrize("tiles", [[0, 1, 1], [1, 1], [True, 2, 2]])
def test_invalid_tiling_is_rejected_before_any_output(tmp_path, tiles):
    folder = folder_fixture(tmp_path)
    config = configuration(tmp_path, [folder], n_tiles=tiles)
    with pytest.raises(ValueError, match="three positive"):
        segment.read_settings(config)
    assert not (folder / "masks").exists()


def test_metadata_json_and_nested_sources_cannot_be_overwritten(tmp_path):
    with pytest.raises(ValueError, match="metadata JSON"):
        segment.read_settings(tmp_path / "._missing.json")
    folder = folder_fixture(tmp_path)
    nested = folder_fixture(folder, "masks")
    config = configuration(tmp_path, [folder, nested])
    with pytest.raises(ValueError, match="contains an input"):
        segment.read_settings(config)
    assert (nested / "processed" / "sample_nuclei.tif").is_file()


def test_model_startup_failure_preserves_previous_results(
    tmp_path, monkeypatch
):
    model = model_fixture(tmp_path, monkeypatch)
    folder = folder_fixture(tmp_path)
    (folder / "masks").mkdir()
    old = folder / "masks" / "old.tif"
    old.write_bytes(b"previous output")
    config = configuration(tmp_path, [folder], model_path=str(model))
    loader = Mock(side_effect=RuntimeError("incompatible weights"))
    monkeypatch.setattr(segment, "load_model", loader)
    monkeypatch.setattr(sys.stdin, "isatty", lambda: False)
    assert segment.main(["-i", str(config)]) == 1
    loader.assert_not_called()
    assert segment.main(["-i", str(config), "--overwrite"]) == 1
    assert old.read_bytes() == b"previous output"


@pytest.mark.parametrize("depth,maximum", [(1, 2), (4, 2), (4, 65536)])
def test_masks_preserve_xyz_calibration_shape_and_large_ids(
    tmp_path, depth, maximum
):
    folder = folder_fixture(tmp_path)
    _, calibration = segment.read_volume(
        folder / "processed" / "sample_nuclei.tif"
    )
    labels = np.zeros((depth, 8, 10), np.uint32)
    labels[:, 2:4, 3:5] = maximum
    labels = segment.compact_dtype(labels)
    path = tmp_path / "mask.tif"
    segment.save_mask(path, labels, calibration)
    with tifffile.TiffFile(path) as tif:
        np.testing.assert_array_equal(tif.asarray(), labels)
        assert tif.series[0].axes == "ZYX"
        metadata = tif.imagej_metadata or tif.shaped_metadata[0]
        assert metadata["spacing"] == 0.75
        assert metadata["unit"] == "um"
        num, den = tif.pages[0].tags["XResolution"].value
        assert den / num == pytest.approx(0.4)
    assert int(tifffile.imread(path).max()) == maximum


def test_constant_and_invalid_volumes_do_not_invent_instances(tmp_path):
    logger = logging.getLogger("test")
    model = Mock()
    labels, _ = segment.predict(np.ones((4, 8, 10)), model, (1, 1, 1), logger)
    assert not labels.any()
    model.predict_instances.assert_not_called()
    data = np.zeros((4, 32, 32))
    data[0, 0, 0] = 1
    with pytest.raises(ValueError, match="percentile range is zero"):
        segment.predict(data, model, (1, 1, 1), logger)
    path = tmp_path / "rgb.tif"
    tifffile.imwrite(path, np.zeros((8, 10, 3), np.uint8), photometric="rgb")
    with pytest.raises(ValueError, match="scalar ZYX"):
        segment.read_volume(path)


def test_batch_failure_continues_without_stale_outputs_or_label_overflow(
    tmp_path, monkeypatch
):
    model_path = model_fixture(tmp_path, monkeypatch)
    folders = [folder_fixture(tmp_path, name) for name in ("first", "second")]
    for folder in folders:
        for name in segment.OUTPUT_DIRS:
            (folder / name).mkdir()
            (folder / name / "old.txt").write_text("stale output")
    bad = folders[0] / "processed" / "broken.tif"
    bad.write_bytes(b"invalid TIFF")
    source = folders[1] / "processed" / "sample_nuclei.tif"
    before = source.read_bytes()
    labels = np.zeros((4, 8, 10), np.uint32)
    labels[1:3, 3:5, 3:5] = 65536
    model = Mock()
    model.predict_instances.return_value = (labels, {})
    tf = Mock()
    loader = Mock(return_value=(model, tf))
    monkeypatch.setattr(segment, "load_model", loader)
    config = configuration(tmp_path, folders, model_path=str(model_path))
    assert segment.main(["-i", str(config), "--overwrite"]) == 1
    loader.assert_called_once()
    tf.keras.backend.clear_session.assert_called_once()
    assert source.read_bytes() == before
    for folder in folders:
        status = json.loads((folder / "masks" / "run_status.json").read_text())
        assert len(status["model_files_sha256"]) == 3
        good = next(r for r in status["images"] if r["status"] == "success")
        assert good["nuclei"] == 1  # ID is 65536, not the number of nuclei.
        assert status["status"] == (
            "failed" if folder == folders[0] else "success"
        )
        assert not list(folder.glob("*/old.txt"))
        assert (folder / "masks" / "sample_nuclei_QC.png").stat().st_size > 0
    assert (
        str(folders[1])
        not in (folders[0] / "nuclei_analysis2.log").read_text()
    )
    assert not (folders[0] / "masks" / "broken_mask.tif").exists()


@pytest.mark.parametrize("option", ["--help", "--version", "--list-models"])
def test_installed_cli_is_lazy_and_models_are_packaged(tmp_path, option):
    script = (
        "import sys; from nuclei_layers_assay.segment import main; "
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
    if option == "--list-models":
        assert len(result.stdout.splitlines()) == 5
        assert "Stardist3D-512-60x" in result.stdout
