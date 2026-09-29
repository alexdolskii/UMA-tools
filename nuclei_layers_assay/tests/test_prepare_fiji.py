"""Native regression tests, including actual process termination."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile

pytestmark = [
    pytest.mark.fiji,
    pytest.mark.skipif(
        os.environ.get("UMA_RUN_FIJI_TESTS") != "1",
        reason="Set UMA_RUN_FIJI_TESTS=1 for native Fiji tests",
    ),
]


def write_stack(path, channels=2):
    """Distinct channels/Z planes catch accidental slice/channel selection."""
    z, c, y, x = np.indices((7, channels, 32, 40))
    data = (200 + 111 * z + 1500 * c + 13 * x + y * y).astype("uint16")
    data[3, -1, 12:20, 12:24] += 8000
    tifffile.imwrite(
        path,
        data,
        imagej=True,
        resolution=(2.5, 1 / 0.75),
        metadata={"axes": "ZCYX", "spacing": 2.0, "unit": "um"},
    )


def command(tmp_path, folders, channel=2):
    config = tmp_path / "input.json"
    config.write_text(
        json.dumps(
            {
                "folder_paths": list(map(str, folders)),
                "nuclei_channel": channel,
                "gaussian_sigma": 4.0,
                "mean_radius": 3,
            }
        )
    )
    executable = Path(sys.executable).with_name("uma_nla_prepare")
    result = subprocess.run(
        [str(executable), "-i", str(config)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=300,
    )
    (tmp_path / "command.log").write_text(result.stdout + result.stderr)
    return result


REFERENCE = r"""
import sys
from pathlib import Path
from uma_tools.imagej import initialize_imagej, shutdown_imagej_workers
from scyjava import jimport

context = initialize_imagej()
source, target, channel = sys.argv[1:]
original = output = None
try:
    IJ = jimport("ij.IJ")
    # Original Bio-Formats defaults for ND2; standard ImageJ TIFF import.
    original = (
        jimport("loci.plugins.BF").openImagePlus(source)[0]
        if source.lower().endswith(".nd2") else IJ.openImage(source)
    )
    macro = (
        'run("Duplicate...", "title=reference duplicate channels=' +
        channel + '");'
        'run("8-bit");'
        'run("Gaussian Blur 3D...", "x=4.0 y=4.0 z=4.0");'
        'run("Mean 3D...", "x=3 y=3 z=3");'
    )
    output = jimport("ij.macro.Interpreter")().runBatchMacro(macro, original)
    assert output is not None, "Reference macro produced no image"
    assert jimport("ij.io.FileSaver")(output).saveAsTiffStack(target)
    cal = original.getCalibration()
    Path(target + ".json").write_text(__import__("json").dumps({
        "dimensions": list(original.getDimensions()),
        "scale": [cal.pixelWidth, cal.pixelHeight, cal.pixelDepth],
        "unit": str(cal.getUnit()),
        "java": str(jimport("java.lang.System").getProperty("java.version")),
    }))
finally:
    for image in (output, original):
        if image is not None:
            image.changes = False
            image.close()
    context.dispose()
    shutdown_imagej_workers()
"""


def reference(tmp_path, source, channel):
    target = tmp_path / "reference.tif"
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            REFERENCE,
            str(source),
            str(target),
            str(channel),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=300,
    )
    (tmp_path / "reference.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    return target


def assert_same_tiff(actual, expected):
    """Compare every voxel and the saved XYZ calibration."""
    with tifffile.TiffFile(actual) as a, tifffile.TiffFile(expected) as b:
        np.testing.assert_array_equal(a.asarray(), b.asarray())
        assert a.asarray().dtype == np.uint8
        assert a.imagej_metadata["spacing"] == b.imagej_metadata["spacing"]
        assert a.imagej_metadata["unit"] == b.imagej_metadata["unit"]
        for name in ("XResolution", "YResolution"):
            assert a.pages[0].tags[name].value == b.pages[0].tags[name].value


def test_two_folders_match_original_filters_and_return(tmp_path):
    folders = [tmp_path / name for name in ("first", "second")]
    for folder in folders:
        folder.mkdir()
        write_stack(folder / "input.tif")
        (folder / "._input.tif").write_bytes(b"metadata, not an image")
    result = command(tmp_path, folders)
    assert result.returncode == 0, result.stdout + result.stderr
    expected = reference(tmp_path, folders[0] / "input.tif", 2)
    for folder in folders:
        outputs = list((folder / "processed").glob("*.tif"))
        assert len(outputs) == 1
        assert tifffile.imread(outputs[0]).shape == (7, 32, 40)
        assert_same_tiff(outputs[0], expected)
        with tifffile.TiffFile(outputs[0]) as tif:
            assert tif.imagej_metadata["spacing"] == 2.0
            xnum, xden = tif.pages[0].tags["XResolution"].value
            ynum, yden = tif.pages[0].tags["YResolution"].value
            assert xden / xnum == pytest.approx(0.4, rel=1e-6)
            assert yden / ynum == pytest.approx(0.75, rel=1e-6)


def test_image_error_returns_failure_after_workers_have_run(tmp_path):
    folder = tmp_path / "source"
    folder.mkdir()
    write_stack(folder / "a_valid.tif")
    write_stack(folder / "z_invalid_channel.tif", channels=1)
    result = command(tmp_path, [folder])
    assert result.returncode == 1, result.stdout + result.stderr
    assert "1 saved; 1 failed" in result.stdout
    assert (folder / "processed" / "a_valid_nuclei.tif").is_file()
    assert not (folder / "processed" / "z_invalid_channel_nuclei.tif").exists()


@pytest.mark.skipif(
    not os.environ.get("UMA_NLA_TEST_ND2"),
    reason="Optional real ND2: set UMA_NLA_TEST_ND2 to its path",
)
def test_real_nd2_channel_one_matches_original_and_returns(tmp_path):
    original = Path(os.environ["UMA_NLA_TEST_ND2"]).resolve()
    folder = tmp_path / "source"
    folder.mkdir()
    source = folder / original.name
    shutil.copy2(original, source)
    result = command(tmp_path, [folder], channel=1)
    assert result.returncode == 0, result.stdout + result.stderr
    expected = reference(tmp_path, source, 1)
    actual = folder / "processed" / f"{source.stem}_nuclei.tif"
    metadata = json.loads(Path(str(expected) + ".json").read_text())
    width, height, _, slices, _ = metadata["dimensions"]
    assert tifffile.imread(actual).shape == (slices, height, width)
    assert_same_tiff(actual, expected)
