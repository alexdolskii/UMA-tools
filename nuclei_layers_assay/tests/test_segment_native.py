"""Opt-in checks with the shipped weights and real TensorFlow processes."""

import json
import os
import shutil
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import tifffile
from nuclei_layers_assay.segment_models import bundled_models

pytestmark = [
    pytest.mark.stardist,
    pytest.mark.skipif(
        os.environ.get("UMA_RUN_STARDIST_TESTS") != "1",
        reason="Set UMA_RUN_STARDIST_TESTS=1 for native model tests",
    ),
]


def process(tmp_path, args, timeout=300):
    env = {
        **os.environ,
        "TF_NUM_INTRAOP_THREADS": "2",
        "TF_NUM_INTEROP_THREADS": "2",
        "OMP_NUM_THREADS": "2",
        "TF_CPP_MIN_LOG_LEVEL": "3",
    }
    result = subprocess.run(
        args,
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    (tmp_path / "last_process.log").write_text(result.stdout + result.stderr)
    return result


def stack(path):
    z, y, x = np.indices((16, 64, 64))
    data = (
        1000
        * np.exp(-((z - 8) ** 2 / 9 + (y - 32) ** 2 / 64 + (x - 32) ** 2 / 64))
    ).astype(np.uint16)
    tifffile.imwrite(path, data, metadata={"axes": "ZYX"})


LOAD_ALL = r"""
import logging
import numpy as np
import tifffile
import sys
from nuclei_layers_assay.segment import load_model, predict
from nuclei_layers_assay.segment_models import bundled_models, validate_model
models = bundled_models()
assert len(models) == 5
volume = tifffile.imread(sys.argv[1])
for name,path in models.items():
    model, tf = load_model(path, validate_model(path))
    labels, _ = predict(volume, model, (1,1,1), logging.getLogger('test'))
    assert labels.shape == volume.shape
    assert labels.dtype == np.uint16
    print('VERIFIED:', name, flush=True)
    tf.keras.backend.clear_session()
"""


def test_all_five_packaged_models_load_predict_and_exit(tmp_path):
    image = tmp_path / "input.tif"
    stack(image)
    result = process(tmp_path, [sys.executable, "-c", LOAD_ALL, str(image)])
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count("VERIFIED:") == 5


def test_cli_continues_after_bad_tiff_and_exits_after_inference(tmp_path):
    folders = [tmp_path / name for name in ("first", "second")]
    for folder in folders:
        (folder / "processed").mkdir(parents=True)
        stack(folder / "processed" / "a_valid_nuclei.tif")
        (folder / "processed" / "._a_valid_nuclei.tif").write_bytes(
            b"metadata"
        )
    (folders[0] / "processed" / "z_broken.tif").write_bytes(b"not a TIFF")
    config = tmp_path / "input.json"
    config.write_text(
        json.dumps(
            {
                "folder_paths": list(map(str, folders)),
                "model_path": "stardist-512-v3-2",
                "n_tiles": [1, 1, 1],
            }
        )
    )
    command = str(Path(sys.executable).with_name("uma_nla_segment"))
    result = process(tmp_path, [command, "-i", str(config)])
    assert result.returncode == 1, result.stdout + result.stderr
    assert "2 saved; 1 failed" in result.stdout
    for index, folder in enumerate(folders):
        masks = folder / "masks"
        assert len(list(masks.glob("*_mask.tif"))) == 1
        status = json.loads((masks / "run_status.json").read_text())
        assert status["status"] == ("failed" if index == 0 else "success")


# Original stage-2 numerical pipeline, kept independent of the new helper.
REFERENCE = r"""
import sys
import numpy as np
import tifffile
from csbdeep.utils import normalize
from stardist.models import StarDist3D
from pathlib import Path
source, target, model_path = map(Path, sys.argv[1:])
np.random.seed(6)
model = StarDist3D(None, name=model_path.name, basedir=str(model_path.parent))
image = normalize(tifffile.imread(source), 1, 99.8, axis=None)
if image.ndim == 2:
    image = image[np.newaxis, ...]
labels, _ = model.predict_instances(image, n_tiles=(2,4,4))
tifffile.imwrite(target, labels.astype(np.uint16))
"""


@pytest.mark.skipif(
    not os.environ.get("UMA_NLA_TEST_PREPARED"),
    reason="Set UMA_NLA_TEST_PREPARED to a representative prepared DAPI TIFF",
)
def test_prepared_dapi_matches_original_instances_and_calibration(tmp_path):
    folder = tmp_path / "source"
    (folder / "processed").mkdir(parents=True)
    source = folder / "processed" / "input_nuclei.tif"
    shutil.copy2(os.environ["UMA_NLA_TEST_PREPARED"], source)
    config = tmp_path / "input.json"
    config.write_text(
        json.dumps(
            {
                "folder_paths": [str(folder)],
                "model_path": "stardist-512-v3-2",
                "n_tiles": [2, 4, 4],
            }
        )
    )
    command = str(Path(sys.executable).with_name("uma_nla_segment"))
    result = process(tmp_path, [command, "-i", str(config)], timeout=900)
    assert result.returncode == 0, result.stdout + result.stderr
    expected = tmp_path / "reference.tif"
    model = bundled_models()["stardist-512-v3-2"]
    result = process(
        tmp_path,
        [
            sys.executable,
            "-c",
            REFERENCE,
            str(source),
            str(expected),
            str(model),
        ],
        timeout=900,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    actual = folder / "masks" / "input_nuclei_mask.tif"
    np.testing.assert_array_equal(
        tifffile.imread(actual), tifffile.imread(expected)
    )
    assert tifffile.imread(actual).max() > 0
    with tifffile.TiffFile(actual) as mask, tifffile.TiffFile(source) as image:
        assert (
            mask.imagej_metadata["spacing"] == image.imagej_metadata["spacing"]
        )
        assert mask.imagej_metadata["unit"] == image.imagej_metadata["unit"]
        for key in ("XResolution", "YResolution"):
            # TIFF may reduce the rational without changing its exact value.
            assert (
                Fraction(*mask.pages[0].tags[key].value)
                == Fraction(*image.pages[0].tags[key].value)
            )
