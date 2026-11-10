"""Export a portable, read-only web view of all recognized UMA runs."""

from __future__ import annotations

import argparse
import json
import platform
import random
import shutil
import sys
import traceback
from pathlib import Path

from . import package_version
from .config import read_config
from .files import assay_directory, safe_label, save_json, sha256_file
from .progress import (
    console,
    folder_logged,
    outcome,
    register_sources,
    update_activity,
)
from .run import close_logger, make_logger, unique_output, utc_now
from .web_report_assets import Assets
from .web_report_data import (
    PREFIX,
    discover,
    link_observation,
    load_run,
    observations,
    truth,
)
from .web_report_plots import charts


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Portable offline report from saved UMA analysis runs"
    )
    parser.add_argument(
        "-i",
        "-I",
        "--input",
        required=True,
        help="JSON containing source folder_paths",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {package_version()}"
    )
    return parser.parse_args(argv)


def write_js(path, target, data):
    """Local classic scripts work under file:// without HTTP requests."""
    payload = json.dumps(
        data, ensure_ascii=True, allow_nan=False, separators=(",", ":")
    )
    path.write_text(f"{target}={payload};\n", encoding="utf-8")


def gallery(run):
    """Stable samples from saved populations; never change their membership."""
    groups = {}
    for index, observation in enumerate(run["observations"]):
        row = observation["row"]
        if run["kind"] != "Report":
            population = "Available intermediate images"
        elif row.get("FIA_Match_Status") not in (None, "", "MATCHED"):
            population = "Unavailable for FIA-matched views"
        elif truth(row.get("Included_In_Filtered_Plots")):
            population = "Retained"
        else:
            population = "Excluded from filtered views"
        groups.setdefault((observation["group"], population), []).append(index)
    rng = random.Random(0)
    result = []
    for (group, population), indices in sorted(groups.items()):
        ordered = sorted(
            indices, key=lambda i: run["observations"][i]["image"]
        )
        sample = sorted(rng.sample(ordered, min(5, len(ordered))))
        result.append(
            {
                "group": group,
                "population": population,
                "indices": sample,
                "total": len(indices),
            }
        )
    return result


def prepare_run(path, kind, output, assets, logger):
    run = load_run(path, kind, output, logger)
    run["charts"] = charts(run)
    run["observations"] = observations(run)
    if run["integrity_ok"]:
        for observation in run["observations"]:
            link_observation(run, observation, assets)
    else:
        run["warnings"].append(
            "Interactive plots and image links disabled: source integrity "
            "could not be verified. Original files remain downloadable."
        )
    run["gallery"] = gallery(run)
    for table in run["tables"].values():
        # The downloadable files retain original decimal precision/encoding.
        table["row_count"] = len(table["rows"])
    # A malformed metadata payload must not abort other saved runs later
    # during JavaScript serialization. Literal CSV strings stay unchanged.
    json.dumps(run, allow_nan=False)
    logger.info(
        "RUN_IMPORTED | %s | %s | %s | images=%d | plots=%d",
        path,
        kind,
        run["status"],
        len(run["observations"]),
        len(run["charts"]),
    )
    return run


def summary(run):
    state = run["documents"].get("run_status.json", {})
    return {
        key: run[key]
        for key in (
            "id",
            "name",
            "kind",
            "status",
            "source",
            "route",
            "warnings",
        )
    } | {
        "images": len(run["observations"]),
        "plots": len(run["charts"]),
        "started": state.get("started_utc", ""),
        "parameters": run["parameters"],
        "script": f"data/{run['id']}.js",
    }


@folder_logged("source")
def export_folder(source, input_json):
    """Keep every prior output; failures in one folder do not stop the next."""
    run_id, output = unique_output(
        assay_directory(source, create=True),
        PREFIX + safe_label(source.name)[:70] + "_",
        include_pid=True,
    )
    logger = make_logger(output)
    state = {
        "schema": "uma-web-report-1",
        "run_id": run_id,
        "status": "RUNNING",
        "started_utc": utc_now(),
        "package_version": package_version(),
        "python": sys.version,
        "platform": platform.platform(),
        "source_folder": str(source),
        "input_json": str(input_json),
        "output": str(output),
        "gallery_seed": 0,
        "gallery_limit_per_population": 5,
        "preview_max_edge": 1024,
        "original_tiffs_copied": False,
        "calculation_policy": (
            "Saved UMA_Report only; no scientific reanalysis"
        ),
    }
    save_json(output / "run_status.json", state)
    logger.info("STARTED | %s", json.dumps(state))
    try:
        sources = discover(source)
        if not sources:
            raise ValueError("No recognizable UMA runs found")
        assets = Assets(source, output, logger)
        runs, failures = [], []
        for index, (path, kind) in enumerate(sources, 1):
            update_activity(f"Loading run {index}/{len(sources)}: {path.name}")
            try:
                runs.append(prepare_run(path, kind, output, assets, logger))
            except Exception as error:
                logger.exception("RUN_IMPORT_FAILED | %s", path)
                failures.append(f"{path.name}: {error}")
        if not runs:
            raise ValueError("No runs could be imported; see run.log")
        previews = assets.build(
            lambda index, total, label: update_activity(
                f"{label}: {index}/{total}"
            )
        )
        data_folder = output / "data"
        data_folder.mkdir()
        catalog = []
        for run in runs:
            for index, observation in enumerate(run["observations"]):
                if observation.get("comparison_key"):
                    catalog.append(
                        {
                            "key": observation["comparison_key"],
                            "run": run["id"],
                            "index": index,
                            "image": observation["image"],
                            "threshold": observation["threshold"],
                            "files": [
                                f
                                for f in observation["files"]
                                if f["label"] == "Area mask"
                            ],
                        }
                    )
            write_js(
                data_folder / f"{run['id']}.js",
                f"window.UMA_RUNS[{json.dumps(run['id'])}]",
                run,
            )
        warnings = failures + assets.failures
        for run in runs:
            warnings.extend(f"{run['name']}: {w}" for w in run["warnings"])
            for observation in run["observations"]:
                warnings.extend(
                    f"{run['name']} / {observation['image']}: {w}"
                    for w in observation["warnings"]
                    if not w.startswith("Unmasked XZ")
                )
        state.update(
            status="PARTIAL" if warnings else "SUCCESS",
            ended_utc=utc_now(),
            recognized_runs=len(sources),
            imported_runs=len(runs),
            warning_count=len(warnings),
            warnings=warnings,
            source_statuses={r["name"]: r["status"] for r in runs},
        )
        payload = {
            "experiment": source.name,
            "state": state,
            "runs": [summary(r) for r in runs],
            "assets": previews,
            "comparisons": catalog,
        }
        write_js(data_folder / "index.js", "window.UMA_INDEX", payload)
        resources = Path(__file__).with_name("web_assets")
        for filename in ("index.html", "app.js", "style.css"):
            shutil.copyfile(resources / filename, output / filename)
        save_json(
            output / "manifest.json",
            {
                **state,
                "display_ranges": assets.ranges,
                "assets": previews,
                "runs": [
                    {
                        "name": r["name"],
                        "source": r["source"],
                        "files_sha256": r["hashes"],
                    }
                    for r in runs
                ],
            },
        )
        save_json(output / "run_status.json", state)
        (output / "README.txt").write_text(
            "Open index.html in a desktop browser. No server is needed.\n"
            "Share this entire folder, not index.html alone.\n"
            "Annotations: save the downloaded HTML beside index.html.\n"
            "TIFF/ND2 files are not copied; their source paths are recorded.\n"
            "Display previews are not measurement images. Saved filtering "
            "and statistics are not recalculated.\n",
            encoding="utf-8",
        )
        integrity = {}
        for file in output.rglob("*"):
            if file.is_file() and file.name != "run.log":
                integrity[file.relative_to(output).as_posix()] = sha256_file(
                    file
                )
        save_json(output / "export_checksums.json", integrity)
        logger.info(
            "%s | %d runs | %d warnings | %s",
            state["status"],
            len(runs),
            len(warnings),
            output / "index.html",
        )
        outcome(
            state["status"],
            f"{len(runs)} runs; {len(warnings)} warnings. "
            f"Open: {output / 'index.html'}",
        )
        return state["status"]
    except (Exception, KeyboardInterrupt) as error:
        logger.exception("WEB_REPORT_FAILED")
        state.update(
            status="CANCELLED"
            if isinstance(error, KeyboardInterrupt)
            else "ERROR",
            error=str(error),
            ended_utc=utc_now(),
        )
        (output / "traceback.txt").write_text(traceback.format_exc())
        save_json(output / "run_status.json", state)
        outcome(state["status"], f"{error}. Diagnostics: {output}")
        if isinstance(error, KeyboardInterrupt):
            raise
        return "ERROR"
    finally:
        close_logger(logger)


def main(argv=None):
    args = parse_args(argv)
    input_json = Path(args.input).expanduser().absolute()
    folders = list(dict.fromkeys(read_config(input_json)))
    register_sources(folders, input_json)
    results = []
    for source in folders:
        try:
            results.append(export_folder(source, input_json))
        except (OSError, ValueError) as error:
            console(f"ERROR: {source}: {error}")
            results.append("ERROR")
    console(
        f"Web reports finished: {results.count('SUCCESS')} succeeded; "
        f"{results.count('PARTIAL')} partial; {results.count('ERROR')} failed."
    )
    return int(any(r != "SUCCESS" for r in results))
