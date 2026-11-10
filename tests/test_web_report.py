"""Offline UMA exports must preserve observations and measurement lineage."""

import json
import logging
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from PIL import Image
from shared_plate_fixture import (
    read_rows,
    seed_all,
    source_fixture,
)

from uma_tools import collect_results, report, web_report
from uma_tools.files import save_csv, sha256_file
from uma_tools.web_report_assets import Assets
from uma_tools.web_report_data import (
    Resolver,
    discover,
    link_nuclei,
    load_run,
    observations,
)
from uma_tools.web_report_plots import charts


def load_script(path):
    return json.loads(path.read_text().split("=", 1)[1].rstrip(";\n"))


def hash_tree(root):
    return {
        str(p.relative_to(root)): sha256_file(p)
        for p in root.rglob("*")
        if p.is_file()
    }


def write_image(path, value=100, mode="signal"):
    path.parent.mkdir(parents=True, exist_ok=True)
    array = np.arange(10000, dtype=np.uint16).reshape(100, 100) % value
    if mode == "binary":
        array = (array > value // 2).astype(np.uint8) * 255
    elif mode == "color":
        array = np.stack([array % 255] * 3, axis=-1).astype(np.uint8)
    Image.fromarray(array).save(path)


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    root = tmp_path_factory.mktemp("web-report")
    source, config, design = source_fixture(root)
    outputs = seed_all(source, design, thresholds=(1500, 2000))
    for output in outputs:
        selection = read_rows(output / "source_selection.csv")
        names = [
            r["File_Name"] for r in selection if r["Status"] == "SELECTED"
        ]
        if output.name.startswith("Alignment"):
            for name in names:
                stem = Path(name).stem + "_processed"
                write_image(output / f"{stem}.tif")
                write_image(
                    output / "Images" / f"{stem}_orientation_composition.png",
                    mode="color",
                )
                write_image(
                    output
                    / "Images"
                    / "normalized_images"
                    / f"{stem}_normalized_orientation.png",
                    mode="color",
                )
        elif output.name.startswith("Thickness"):
            for name in names:
                write_image(output / f"Mask_{name}.tif", mode="binary")
                write_image(output / f"Local_Thickness_{name}.tif")
        else:
            path = output / "Fibronectin_Area_Summary.csv"
            records = read_rows(path)
            for row in records:
                name = row["File_Name"]
                projection = (
                    output / "Projections_32bit" / f"{name}_FN_SUM32.tif"
                )
                mask = output / "Masks" / f"{name}_mask.tif"
                write_image(projection, 200)
                write_image(mask, int(row["Threshold_Lower"]), "binary")
                row.update(
                    Source_Projection_Path=str(projection),
                    Projection_SHA256=sha256_file(projection),
                    Mask_Path=str(mask),
                    Mask_SHA256=sha256_file(mask),
                    Source_Original_Path=str(source / name),
                    Source_SHA256=sha256_file(source / name),
                )
            save_csv(path, list(records[0]), records)
            state = json.loads((output / "run_status.json").read_text())
            # Match ImageRun's completion record after adding saved artifacts.
            state["summary_sha256"] = sha256_file(path)
            (output / "run_status.json").write_text(json.dumps(state))
    assert collect_results.main(["-i", str(config), "--non-interactive"]) == 0
    assert (
        report.main(
            [
                "-i",
                str(config),
                "--non-interactive",
                "--stats-unit",
                "well",
                "--plot-format",
                "both",
            ]
        )
        == 0
    )
    before = hash_tree(source)
    assert web_report.main(["-i", str(config)]) == 0
    output = next((source / "uma_assay").glob("UMA_Web_Report_*"))
    assert all(
        sha256_file(source / relative) == digest
        for relative, digest in before.items()
    )
    index = load_script(output / "data" / "index.js")
    runs = [load_script(output / r["script"]) for r in index["runs"]]
    return source, config, output, index, runs


def test_complete_pipeline_population_boxes_statistics_and_colors(exported):
    source, _, output, index, runs = exported
    assert index["state"]["status"] == "SUCCESS"
    assert len([r for r in runs if r["kind"] == "Report"]) == 6
    for run in runs:
        if run["kind"] != "Report":
            assert not run["charts"]
            continue
        assert len(run["charts"]) == 6
        original = json.loads(
            (Path(run["source"]) / "plot_manifest.json").read_text()
        )
        for plot, saved in zip(run["charts"], original):
            assert {k: v for k, v in plot.items() if k != "points"} == saved
            assert {p["Image_ID"] for p in plot["points"]} == set(
                saved["plotted_image_ids"]
            )
            if "No images" in saved["group_counts"]:
                assert saved["group_counts"]["No images"] == 0
        assert any(c["statistical_comparisons"] for c in run["charts"])
        for record in run["downloads"]:
            assert sha256_file(output / record["path"]) == record["sha256"]
        assert {o["image"] for o in run["observations"]} == {
            r["Image_ID"]
            for r in read_rows(Path(run["source"]) / "merged_data.csv")
        }
    assert not list(output.rglob("*.tif"))
    assert not list(output.rglob("*.nd2"))
    assert len(discover(source)) == len(runs)


def test_thresholds_channels_and_exact_run_images_stay_separate(exported):
    _, _, _, index, runs = exported
    reports = [r for r in runs if r["kind"] == "Report"]
    for run in reports:
        selected = run["documents"]["Inputs/collector_run_status.json"][
            "selected"
        ]
        for observation in run["observations"]:
            for file in observation["files"]:
                asset = index["assets"][file["asset"]]
                assert asset["status"] == "AVAILABLE"
                kind = (
                    "Area"
                    if file["label"].startswith("Area")
                    else (
                        "Thickness"
                        if "thickness" in file["label"].lower()
                        else "Alignment"
                    )
                )
                assert Path(asset["source"]).is_relative_to(
                    selected[kind]["run_directory"]
                )
            candidates = [
                c
                for c in index["comparisons"]
                if c["key"] == observation["comparison_key"]
            ]
            assert {c["threshold"]["Threshold_Lower"] for c in candidates} == {
                "1500",
                "2000",
            }
            for candidate in candidates:
                other = next(r for r in runs if r["id"] == candidate["run"])
                assert other["route"] == run["route"]


def test_gallery_is_reproducible_and_keeps_excluded_population(exported):
    *_, runs = exported
    for run in runs:
        assert run["gallery"] == web_report.gallery(run)
        assert all(len(g["indices"]) <= 5 for g in run["gallery"])
        if run["kind"] == "Report":
            assert "Excluded from filtered views" in {
                g["population"] for g in run["gallery"]
            }


def test_unsafe_missing_and_relocated_paths(tmp_path):
    source = tmp_path / "source"
    run = source / "uma_assay" / "Area_assay_results_1"
    run.mkdir(parents=True)
    target = run / "mask.tif"
    target.write_bytes(b"image")
    resolver = Resolver(source)
    assert (
        resolver.resolve("C:\\old\\uma_assay\\Area_assay_results_1\\mask.tif")
        == target
    )
    assert (
        resolver.resolve("/old/uma_assay/Area_assay_results_1/mask.tif")
        == target
    )
    assert (
        resolver.resolve("../source/uma_assay/Area_assay_results_1/mask.tif")
        is None
    )
    other = tmp_path / "foreign.tif"
    other.write_bytes(b"foreign")
    (run / "linked.tif").symlink_to(other)
    assert resolver.resolve("linked.tif", run) is None
    assert resolver.resolve("._mask.tif", run) is None
    assert resolver.resolve(str(other)) is None
    assert resolver.resolve("/unknown/mask.tif") is None


def test_preview_common_ranges_stack_rejection_and_hash_checks(tmp_path):
    run = tmp_path / "uma_assay" / "Area_assay_results_1"
    run.mkdir(parents=True)
    assets = Assets(tmp_path, tmp_path / "export", logging.getLogger("assets"))
    for name, high in (("a", 100), ("b", 200)):
        write_image(run / f"{name}.tif", high)
        assets.add(f"{name}.tif", run, name, "same marker")
    image = Image.fromarray(np.ones((5, 5), np.uint8))
    image.save(run / "stack.tif", save_all=True, append_images=[image])
    assets.add("stack.tif", run, "stack", "same marker")
    assets.add("a.tif", run, "wrong hash", "same marker", expected="0" * 64)
    assets.add("missing.tif", run, "missing", "same marker")
    result = assets.build(lambda *_: None)
    available = [a for a in result.values() if a["status"] == "AVAILABLE"]
    assert len(available) == 2
    assert all(a["display_range"] == [0, 199] for a in available)
    assert len(assets.failures) == 3


def test_modified_report_does_not_generate_plots(exported, tmp_path):
    _, _, _, _, runs = exported
    original = next(r for r in runs if r["kind"] == "Report")
    path = tmp_path / "source" / "uma_assay" / "UMA_Report_tampered"
    shutil.copytree(original["source"], path)
    table = path / "plot_data.csv"
    table.write_text(table.read_text() + "unexpected,row\n")
    loaded = load_run(
        path, "Report", tmp_path / "export", logging.getLogger("load")
    )
    assert not loaded["integrity_ok"]
    assert any("checksum" in w for w in loaded["warnings"])
    assert charts(loaded) == []
    assert observations(loaded)


def test_fia_uses_archived_run_and_requires_projection_checksum(tmp_path):
    parent = tmp_path / "fia_assay" / "Final_Nuclei_Mask_1"
    parent.mkdir(parents=True)
    meta = {"stardist_folder": "not available"}
    (parent / "nuclei_run.json").write_text(json.dumps(meta))
    run = {
        "documents": {
            "fia_selection.json": {"run": str(parent)},
            "Inputs/FIA/01_nuclei_run.json": meta,
        }
    }
    row = {
        "FIA_Match_Status": "MATCHED",
        "FIA_Run_ID": parent.name,
        "FIA_Mask_Name": (
            "one_nuclei_projection_StarDist_processed_processed.tif"
        ),
    }
    observation = {"row": row, "warnings": []}
    calls = []
    assets = Assets(tmp_path, tmp_path / "export", logging.getLogger("fia"))
    link_nuclei(
        run, observation, assets, lambda *args, **kw: calls.append(args)
    )
    assert len(calls) == 1 and calls[0][2] == "FIA nuclear mask"
    assert "source hash" in observation["warnings"][0]
    (parent / "nuclei_run.json").write_text("{}")
    calls.clear()
    link_nuclei(
        run, observation, assets, lambda *args, **kw: calls.append(args)
    )
    assert not calls


def test_fia_projection_requires_exact_recorded_bytes(tmp_path):
    parent = tmp_path / "fia_assay" / "Final_Nuclei_Mask_1"
    star = tmp_path / "fia_assay" / "StarDist_1"
    parent.mkdir(parents=True)
    star.mkdir()
    meta = {"stardist_folder": str(star)}
    (parent / "nuclei_run.json").write_text(json.dumps(meta))
    name = "one_nuclei_projection.tif"
    projection = tmp_path / "fia_assay" / "Nuclei" / name
    write_image(projection)
    (star / "stardist_run.json").write_text(
        json.dumps({"sources": {name: sha256_file(projection)}})
    )
    mask = "one_nuclei_projection_StarDist_processed_processed.tif"
    write_image(parent / mask, mode="binary")
    run = {
        "documents": {
            "fia_selection.json": {"run": str(parent)},
            "Inputs/FIA/01_nuclei_run.json": meta,
        }
    }
    observation = {
        "row": {
            "FIA_Match_Status": "MATCHED",
            "FIA_Run_ID": parent.name,
            "FIA_Mask_Name": mask,
        },
        "warnings": [],
    }
    assets = Assets(tmp_path, tmp_path / "export", logging.getLogger("fia"))

    def add(value, base, label, mode="linear", expected=None):
        return assets.add(value, base, label, "fia", mode, expected)

    link_nuclei(run, observation, assets, add)
    assert not observation["warnings"]
    write_image(projection, 200)
    result = assets.build(lambda *_: None)
    assert (
        next(
            a
            for a in result.values()
            if a["label"] == ("FIA nuclei projection")
        )["status"]
        == "UNAVAILABLE"
    )
    assert (
        next(a for a in result.values() if a["label"] == ("FIA nuclear mask"))[
            "status"
        ]
        == "AVAILABLE"
    )


def test_three_hundred_images_and_two_thresholds(tmp_path):
    source = tmp_path / "300 images"
    source.mkdir()
    config = tmp_path / "paths.json"
    config.write_text(json.dumps({"folder_paths": [str(source)]}))
    for lower in (1500, 2000):
        run = source / "uma_assay" / f"Area_assay_results_{lower}"
        run.mkdir(parents=True)
        records = []
        for index in range(300):
            image = f"source_{700 - index}.nd2"
            path = run / f"{image}_projection.tif"
            write_image(path)
            mask = run / f"{image}_mask.tif"
            write_image(mask, mode="binary")
            records.append(
                {
                    "Image_ID": image,
                    "Source_Projection_Path": str(path),
                    "Mask_Path": str(mask),
                    "Threshold_Lower": lower,
                    "Source_SHA256": str(index).zfill(64),
                }
            )
        save_csv(
            run / "Fibronectin_Area_Summary.csv", list(records[0]), records
        )
        (run / "run_status.json").write_text(
            json.dumps(
                {
                    "status": "PARTIAL" if lower == 1500 else "SUCCESS",
                    "fiber_route": {"marker_id": "FN", "source_channel": 3},
                }
            )
        )
    assert web_report.main(["-i", str(config)]) == 0
    output = next((source / "uma_assay").glob("UMA_Web_Report_*"))
    data = load_script(output / "data" / "index.js")
    assert len(data["assets"]) == 1200
    assert len(data["comparisons"]) == 600
    assert all(r["images"] == 300 and r["plots"] == 0 for r in data["runs"])
    assert all(a["status"] == "AVAILABLE" for a in data["assets"].values())


def test_missing_preview_is_partial_and_later_folder_continues(tmp_path):
    folders = []
    for suffix in ("one", "two"):
        source = tmp_path / suffix
        run = source / "uma_assay" / "Area_assay_results_1"
        run.mkdir(parents=True)
        save_csv(
            run / "Fibronectin_Area_Summary.csv",
            ["Image_ID", "Mask_Path"],
            [{"Image_ID": "image.nd2", "Mask_Path": str(run / "missing.tif")}],
        )
        folders.append(str(source))
    config = tmp_path / "paths.json"
    config.write_text(json.dumps({"folder_paths": folders}))
    assert web_report.main(["-i", str(config)]) == 1
    for folder in folders:
        output = next((Path(folder) / "uma_assay").glob("UMA_Web_Report_*"))
        assert (output / "index.html").is_file()
        assert (
            json.loads((output / "run_status.json").read_text())["status"]
            == "PARTIAL"
        )


def test_metadata_json_is_rejected_before_opening(tmp_path):
    with pytest.raises(ValueError, match="metadata"):
        web_report.main(["-i", str(tmp_path / "._config.json")])


@pytest.mark.skipif(
    os.environ.get("UMA_RUN_WEB_UI_TESTS") != "1",
    reason="Set UMA_RUN_WEB_UI_TESTS=1 with Node/jsdom installed",
)
def test_offline_dom_and_portable_folder(exported, tmp_path):
    _, _, output, _, _ = exported
    relocated = tmp_path / "renamed portable report"
    shutil.copytree(output, relocated)
    result = subprocess.run(
        [
            "node",
            str(Path(__file__).with_name("web_report_ui.cjs")),
            str(relocated),
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
