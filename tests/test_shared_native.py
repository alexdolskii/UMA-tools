"""Real Fiji: scheme routing preserves the existing pixels and measurements."""

import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(
    os.environ.get("UMA_RUN_IMAGEJ_TESTS") != "1",
    reason="Opt-in native Fiji check",
)
def test_identical_fiber_pixels_in_c2_and_c3_match_unrouted_engines(tmp_path):
    result = subprocess.run(
        [sys.executable, __file__, str(tmp_path)],
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def _verify_shared_pixels(tmp_path):
    import numpy as np
    import openpyxl
    import tifffile
    from shared_plate_fixture import read_rows, scheme

    from uma_tools import alignment_analysis as alignment
    from uma_tools import area_analysis as area
    from uma_tools import thickness_analysis as thickness
    from uma_tools.area_imagej import threshold_settings
    from uma_tools.fiber_routes import routes, use_route
    from uma_tools.imagej import initialize_imagej, shutdown_imagej_workers
    from uma_tools.plate_design import read_design

    source, reference = tmp_path / "scheme", tmp_path / "unrouted"
    source.mkdir()
    reference.mkdir()
    design = scheme(
        source / "FIA_plate.xlsx",
        channels={"B02": ("Control", 2), "B03": ("Treatment", 3)},
        layout=["Treatment", "Control"],
    )
    book = openpyxl.load_workbook(design.path)
    for col in (4, 5):
        book["Channels"].cell(13, col, "Ignore")
    book.save(design.path)
    book.close()
    design = read_design(design.path)
    fibers = np.zeros((7, 32, 32), dtype=np.uint16)
    fibers[1:6, 5:27, 8:12] = 30000
    fibers[2:5, 9:23, 21:25] = 45000
    for well, channel in [("B02", 2), ("B03", 3)]:
        image = np.zeros((7, 4, 32, 32), dtype=np.uint16)
        image[:, channel - 1] = fibers
        name = f"fiber_Well{well}.tiff"
        tifffile.imwrite(
            source / name,
            image,
            imagej=True,
            resolution=(2, 2),
            metadata={"axes": "ZCYX", "spacing": 0.75, "unit": "um"},
        )
        if channel == 2:
            (reference / name).write_bytes((source / name).read_bytes())
    ij = initialize_imagej()
    holder, outputs = [None], []
    try:
        classes = thickness.import_java_classes()

        def measure(folder, channel):
            assert (
                alignment.process_folder(
                    str(folder), channel, 15, 500, 500, ij
                )
                == "SUCCESS"
            )
            IJ, _, windows, table, duplicator, _ = classes
            assert (
                thickness.process_single_folder(
                    ij,
                    IJ,
                    windows,
                    duplicator,
                    table,
                    str(folder),
                    None,
                    channel,
                )
                == "SUCCESS"
            )
            assert (
                area.process_folder(
                    folder,
                    channel,
                    "sum",
                    threshold_settings(2000, None),
                    holder,
                    None,
                    outputs,
                ).value
                == "SUCCESS"
            )

        measure(reference, 2)
        for route in routes(source, design):
            with use_route(route):
                measure(source, route.channel)
        fields = {
            "Alignment": (
                "Analysis/Alignment_Summary.csv",
                ["Percentage_Fibers_Aligned_Within_15_Degree"],
            ),
            "Thickness": (
                "Thickness_Summary.csv",
                ["Area", "StdDev", "Min", "Max", "Median"],
            ),
            "Area": (
                "Fibronectin_Area_Summary.csv",
                ["FN_Positive_Pixels", "FN_Area_Percent", "FN_Area"],
            ),
        }
        for prefix, (summary, columns) in fields.items():
            expected_folder = next(
                (reference / "uma_assay").glob(prefix + "_assay_results_*")
            )
            expected = read_rows(expected_folder / summary)[0]
            actual_folders = list(
                (source / "uma_assay").glob(prefix + "_assay_results_*")
            )
            assert len(actual_folders) == 2
            for folder in actual_folders:
                rows = read_rows(folder / summary)
                assert len(rows) == 1
                assert {field: float(rows[0][field]) for field in columns} == {
                    field: float(expected[field]) for field in columns
                }
                assert int(rows[0]["Source_Channel"]) in (2, 3)
                for pattern in (
                    "Mask_*.tif",
                    "Local_Thickness_*.tif",
                    "Masks/*.tif",
                    "Projections_32bit/*.tif",
                ):
                    expected_images = list(expected_folder.glob(pattern))
                    if expected_images:
                        actual_image = next(folder.glob(pattern))
                        np.testing.assert_array_equal(
                            tifffile.imread(actual_image),
                            tifffile.imread(expected_images[0]),
                        )
    finally:
        if holder[0] is not None:
            holder[0].close()
        ij.dispose()
        shutdown_imagej_workers()


if __name__ == "__main__":
    _verify_shared_pixels(Path(sys.argv[1]))
