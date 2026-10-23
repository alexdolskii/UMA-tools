"""Read FIA final-nuclei tables and join them to UMA source images.

This is a file-format adapter. It never imports FIA, starts Fiji, reads
marker intensities, or combines nuclei measured by different runs.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from xml.etree.ElementTree import ParseError
from zipfile import BadZipFile

from .report_inputs import save_csv, save_json
from .report_schema import (
    FN_INCLUDED_FLAG,
    FN_LOW_FLAG,
    FN_REASON_COLUMN,
    ValidationError,
)
from .report_validation import filter_counts

COUNT = "Non_border_nuclei_count"
RATIO = "FN_Area_Per_Nucleus (µm²/nucleus)"
LOW = "Below_Nuclei_Threshold"
RUN_PATTERN = re.compile(r"Final_Nuclei_Mask_(\d{8}_\d{6})")
MASK_SUFFIX = "_nuclei_projection_StarDist_processed_processed"
COUNTS = ("Nuclei_count_total", "Border_nuclei_count", COUNT)
SCALES = {
    "um": 1.0,
    "micron": 1.0,
    "microns": 1.0,
    "micrometer": 1.0,
    "micrometers": 1.0,
    "micrometre": 1.0,
    "micrometres": 1.0,
    "nm": 0.001,
    "mm": 1000.0,
    "cm": 10000.0,
    "m": 1000000.0,
}
ADDED_COLUMNS = [
    "FIA_Run_ID",
    "FIA_Mask_Name",
    *COUNTS,
    "Min_Nuclei",
    LOW,
    "FIA_Match_Status",
    "FIA_Identity_Check",
    "FN_Area_um2",
    RATIO,
    "FN_Per_Nucleus_Reason",
]
MATCH_COLUMNS = ["Image_ID", "FIA_Mask_Name", "Status", "Reason"]


def _number(value, label, integer=False):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"Invalid {label}: {value!r}") from error
    if (
        not math.isfinite(result)
        or result < 0
        or (integer and not result.is_integer())
    ):
        raise ValidationError(f"Invalid {label}: {value!r}")
    return int(result) if integer else result


class _Reader:
    """Keep immutable source bytes for validation and later archiving."""

    def __init__(self):
        self.files = {}

    def read(self, path):
        if (
            path.is_symlink()
            or not path.is_file()
            or path.name.startswith(".")
        ):
            raise ValidationError(f"Missing regular FIA input: {path}")
        if path not in self.files:
            self.files[path] = path.read_bytes()
        return self.files[path]

    def json(self, path):
        value = json.loads(self.read(path).decode("utf-8-sig"))
        if not isinstance(value, dict):
            raise ValidationError(f"Expected a JSON object: {path}")
        return value

    def table(self, path, required):
        reader = csv.reader(
            io.StringIO(self.read(path).decode("utf-8-sig")), strict=True
        )
        columns = next(reader, [])
        if (
            not columns
            or any(not c for c in columns)
            or len(set(columns)) != len(columns)
            or not set(required) <= set(columns)
        ):
            raise ValidationError(f"Missing or duplicate FIA columns: {path}")
        rows = []
        for index, values in enumerate(reader, 2):
            if len(values) != len(columns):
                raise ValidationError(f"Invalid CSV row {index}: {path}")
            rows.append(dict(zip(columns, values)))
        return columns, rows


def _verify_excel(reader, path, tables):
    from openpyxl import load_workbook

    book = load_workbook(
        io.BytesIO(reader.read(path)), read_only=True, data_only=False
    )
    try:
        for name, (columns, rows) in tables.items():
            if name not in book:
                raise ValidationError(f"Missing FIA worksheet: {path}/{name}")
            values = list(book[name].values)
            if (
                not values
                or list(values[0]) != columns
                or len(values) != len(rows) + 1
            ):
                raise ValidationError(
                    f"FIA CSV/workbook structure mismatch: {name}"
                )
            for expected, actual in zip(rows, values[1:]):
                if len(actual) > len(columns):
                    raise ValidationError(f"Extra FIA workbook cells: {name}")
                actual = tuple(actual) + (None,) * (len(columns) - len(actual))
                for column, value in zip(columns, actual):
                    text = expected[column]
                    if value is None and text == "":
                        continue
                    if str(value) == text:
                        continue
                    if (
                        name == "Run_Info"
                        and column == "Value"
                        and len(text) > 32767
                        and str(value) == text[:32767]
                    ):
                        # FIA keeps full CSV text past Excel's cell limit.
                        continue
                    try:
                        equal = math.isclose(
                            float(value),
                            float(text),
                            rel_tol=1e-12,
                            abs_tol=1e-12,
                        )
                    except (ValueError, TypeError):
                        equal = False
                    if not equal:
                        raise ValidationError(
                            f"FIA CSV/workbook value mismatch: {name}/{column}"
                        )
    finally:
        book.close()


def _exclusions(reader, run):
    records = {}
    for path in (
        run / "excluded_images.json",
        run.parent / "excluded_images.json",
    ):
        if not path.exists():
            continue
        value = reader.json(path)
        if value.get("schema_version") != 1 or not isinstance(
            value.get("images"), dict
        ):
            raise ValidationError(f"Invalid FIA exclusion registry: {path}")
        for key, row in value["images"].items():
            if (
                not isinstance(row, dict)
                or row.get("Image_key") != key
                or not row.get("Reason")
            ):
                raise ValidationError(f"Invalid FIA exclusion: {path}: {key}")
            records[key] = row
    return records


def load_run(run, source):
    """Validate the shared nuclei-only/full-marker source format."""
    reader = _Reader()
    from .fiber_routes import active_route

    route = active_route()
    if route is not None:
        from .plate_design import read_design

        scheme_meta = reader.json(run / "plate_design.json")
        content = reader.read(run / "plate_design.xlsx")
        scheme = read_design(run / "plate_design.xlsx")
        if (
            hashlib.sha256(content).hexdigest()
            != scheme_meta.get("Workbook_SHA256")
            or scheme_meta.get("Schema_version") != 1
            or scheme_meta.get("Wells") != scheme.wells
            or scheme_meta.get("Panels") != scheme.panels
            or scheme.scope_hash("Nuclei") != scheme_meta.get("Nuclei_SHA256")
            or scheme.scope_hash("Nuclei") != route.design.scope_hash("Nuclei")
        ):
            raise ValidationError(
                "FIA nuclei scheme is missing, changed, or incompatible with "
                "current Nuclei assignments"
            )
    meta = reader.json(run / "nuclei_run.json")
    if (
        meta.get("schema_version") != 1
        or meta.get("status") != "complete"
        or meta.get("morphology_status") != "complete"
        or meta.get("skipped_files")
    ):
        raise ValidationError("Nuclei and morphology must both be complete")
    identity = [
        "Dataset_path",
        "Run_ID",
        "Mask_name",
        "Image_name",
        "Particle_size_px2",
    ]
    tables = {
        "Images": reader.table(
            run / "Nuclei_Images.csv", identity + ["Status", "Error", *COUNTS]
        ),
        "Nuclei": reader.table(
            run / "Nuclei_Morphology.csv",
            identity + ["Nucleus_ID", "Touches_border"],
        ),
        "Run_Info": reader.table(
            run / "Nuclei_Run_Info.csv", ["Parameter", "Value"]
        ),
    }
    info_rows = tables["Run_Info"][1]
    info = {row["Parameter"]: row["Value"] for row in info_rows}
    particle = _number(
        meta.get("particle_size_pixels_squared"), "particle size"
    )
    if (
        len(info) != len(info_rows)
        or info.get("Status") != "complete"
        or info.get("Run_ID") != run.name
        or _number(info.get("Particle_size_px2"), "particle size") != particle
    ):
        raise ValidationError("Inconsistent FIA Run_Info")
    images, nucleus_ids = {}, set()
    populations = Counter()
    for sheet in ("Images", "Nuclei"):
        for row in tables[sheet][1]:
            mask = row["Mask_name"]
            if (
                Path(mask).name != mask
                or "\\" in mask
                or mask.startswith(".")
                or not mask.endswith(MASK_SUFFIX + ".tif")
            ):
                raise ValidationError(
                    f"Unrecognized FIA mask identity: {mask}"
                )
            expected_image = (
                mask.removesuffix("_StarDist_processed_processed.tif") + ".tif"
            )
            if row["Image_name"] != expected_image:
                raise ValidationError(
                    f"FIA image/mask identity mismatch: {mask}"
                )
            if (
                row["Run_ID"] != run.name
                or Path(row["Dataset_path"]).resolve() != source.resolve()
                or _number(row["Particle_size_px2"], "particle size")
                != particle
            ):
                raise ValidationError(
                    f"FIA dataset/run/parameters mismatch: {mask}"
                )
            if sheet == "Images":
                if (
                    mask in images
                    or row["Status"] != "complete"
                    or row["Error"].strip()
                ):
                    raise ValidationError(
                        f"Duplicate or unsuccessful FIA image: {mask}"
                    )
                counts = {key: _number(row[key], key, True) for key in COUNTS}
                if counts[COUNTS[0]] - counts[COUNTS[1]] != counts[COUNT]:
                    raise ValidationError(f"Inconsistent FIA counts: {mask}")
                images[mask] = {**row, **counts}
            else:
                key = (mask, _number(row["Nucleus_ID"], "Nucleus_ID", True))
                if (
                    key in nucleus_ids
                    or key[1] < 1
                    or row["Touches_border"].lower() not in ("false", "0")
                ):
                    raise ValidationError(
                        f"Invalid or duplicate non-border nucleus: {key}"
                    )
                nucleus_ids.add(key)
                populations[mask] += 1
    if (
        not images
        or any(populations[mask] != row[COUNT] for mask, row in images.items())
        or set(populations) - set(images)
    ):
        raise ValidationError(
            "FIA nucleus rows do not reconcile with image counts"
        )
    processed = meta.get("processed_files")
    if (
        not isinstance(processed, list)
        or not all(isinstance(name, str) for name in processed)
        or len(processed) != len(images)
        or {Path(name).stem + "_processed.tif" for name in processed}
        != set(images)
    ):
        raise ValidationError(
            "FIA processed-files manifest differs from image tables"
        )
    _verify_excel(reader, run / "Nuclei_Morphology.xlsx", tables)
    calibration = {}
    if (run / "spatial_calibration.json").exists():
        calibration = reader.json(run / "spatial_calibration.json")
        for mask, snapshot in calibration.items():
            if not isinstance(snapshot, dict):
                raise ValidationError(
                    f"Invalid FIA calibration record: {mask}"
                )
            shape = snapshot.get("Shape_YX")
            if (
                not isinstance(shape, list)
                or len(shape) != 2
                or any(_number(v, "dimension", True) < 1 for v in shape)
            ):
                raise ValidationError(f"Invalid FIA dimensions: {mask}")
    for mask, image in images.items():
        if image.get("Calibration_status") == "calibrated":
            for axis in ("X", "Y"):
                if (
                    _number(image.get(f"Pixel_size_{axis}_um"), "pixel size")
                    <= 0
                ):
                    raise ValidationError(f"Invalid FIA pixel size: {mask}")
    return {
        "run": run,
        "images": images,
        "files": reader.files,
        "particle_size": particle,
        "exclusions": _exclusions(reader, run),
        "calibration": calibration,
    }


def select_run(source, log):
    assay = source / "fia_assay"
    if not assay.is_dir() or assay.is_symlink():
        raise ValidationError(f"No regular FIA results directory: {assay}")
    for run in sorted(assay.iterdir(), reverse=True):
        match = RUN_PATTERN.fullmatch(run.name)
        if not match or not run.is_dir() or run.is_symlink():
            continue
        try:
            datetime.strptime(match[1], "%Y%m%d_%H%M%S")
            bundle = load_run(run, source)
        except (
            OSError,
            ValueError,
            csv.Error,
            BadZipFile,
            ParseError,
            ValidationError,
        ) as error:
            log.event("WARNING", "FIA source rejected", f"{run}: {error}")
            continue
        log.event(
            "INFO",
            "FIA source selected",
            f"{run}; particle size {bundle['particle_size']:g} px²",
        )
        return bundle
    raise ValidationError(f"No completed valid FIA nuclei run in {assay}")


def _archive(bundle, directory, manifests):
    folder = directory / "Inputs" / "FIA"
    folder.mkdir()
    for index, (source, content) in enumerate(bundle["files"].items(), 1):
        target = folder / f"{index:02d}_{source.name}"
        target.write_bytes(content)
        manifests.append(
            {
                "Input": "fia_" + str(index),
                "Path": str(source),
                "Archived_Path": str(target),
                "Archive_Relative_Path": str(target.relative_to(directory)),
                "Bytes": len(content),
                "SHA256": hashlib.sha256(content).hexdigest(),
            }
        )
    save_csv(directory / "input_manifest.csv", list(manifests[0]), manifests)


def _unit(value):
    return str(value).strip().lower().replace("µ", "u").replace("μ", "u")


def _identity(row, image, snapshot):
    """Reject contradictory geometry; expose missing evidence explicitly."""
    checks = ["exact source stem", "dataset", "run"]
    width, height = (
        row.get("FN_Source__" + key)
        for key in ("Width_Pixels", "Height_Pixels")
    )
    if snapshot.get("Shape_YX") and width is not None and height is not None:
        expected = [
            _number(height, "height", True),
            _number(width, "width", True),
        ]
        if expected != snapshot["Shape_YX"]:
            raise ValidationError("FIA/UMA image dimensions differ")
        checks.append("XY dimensions")
    else:
        checks.append("XY dimensions not recorded in both inputs")
    scale = SCALES.get(_unit(row.get("FN_Source__Pixel_Unit")))
    if scale and image.get("Calibration_status") == "calibrated":
        for axis, field in (("X", "Pixel_Width"), ("Y", "Pixel_Height")):
            uma = _number(row.get("FN_Source__" + field), field) * scale
            fia = _number(image.get(f"Pixel_size_{axis}_um"), "FIA pixel size")
            if (
                uma <= 0
                or fia <= 0
                or not math.isclose(uma, fia, rel_tol=1e-5, abs_tol=1e-9)
            ):
                raise ValidationError("FIA/UMA physical pixel sizes differ")
        checks.append("XY calibration")
    else:
        checks.append("XY calibration not recorded in both inputs")
    original_hash = image.get("Source_SHA256") or snapshot.get(
        "Original_SHA256"
    )
    uma_hash = row.get("FN_Source__Source_SHA256")
    if original_hash and uma_hash:
        if original_hash.lower() != str(uma_hash).lower():
            raise ValidationError("FIA/UMA original file hashes differ")
        checks.append("original SHA256")
    else:
        checks.append("original SHA256 not recorded in both inputs")
    return "; ".join(checks)


def _area_um2(row):
    area = row.get("FN_Source__FN_Area")
    unit = _unit(row.get("FN_Source__Area_Unit")).replace("²", "^2")
    scale = (
        SCALES.get(unit.removesuffix("^2")) if unit.endswith("^2") else None
    )
    if area in (None, "") or scale is None:
        return None
    value = _number(area, "FN_Area") * scale**2
    pixel_unit = SCALES.get(_unit(row.get("FN_Source__Pixel_Unit")))
    if pixel_unit and all(
        row.get("FN_Source__" + key) is not None
        for key in ("FN_Positive_Pixels", "Pixel_Width", "Pixel_Height")
    ):
        positive, width, height = [
            _number(row["FN_Source__" + key], key)
            for key in ("FN_Positive_Pixels", "Pixel_Width", "Pixel_Height")
        ]
        if (
            width <= 0
            or height <= 0
            or not math.isclose(
                value,
                positive * width * height * pixel_unit**2,
                rel_tol=1e-6,
                abs_tol=1e-8,
            )
        ):
            raise ValidationError(
                "FN physical area disagrees with positive pixels "
                "and calibration"
            )
    return value


def attach_nuclei(
    data, source, directory, manifests, minimum, log, *, bundle=None
):
    """Attach counts and audit every image before filtering matched rows."""
    bundle = select_run(source, log) if bundle is None else bundle
    _archive(bundle, directory, manifests)
    lookup = {}
    for mask, row in bundle["images"].items():
        stem = mask.removesuffix(MASK_SUFFIX + ".tif")
        if stem in lookup:
            raise ValidationError(f"Ambiguous FIA original stem: {stem}")
        lookup[stem] = (mask, row)
    seen, matches, audit, accepted = set(), [], [], []
    for row in data["rows"]:
        stem = Path(row["Image_ID"]).stem
        if stem in seen:
            raise ValidationError(f"Ambiguous UMA original stem: {stem}")
        seen.add(stem)
        row.update(dict.fromkeys(ADDED_COLUMNS))
        row.update(FIA_Run_ID=bundle["run"].name, Min_Nuclei=minimum)
        item = lookup.get(stem)
        status, reason = "MATCHED", ""
        if stem in bundle["exclusions"]:
            status, reason = (
                "PROCESSING_EXCLUDED",
                bundle["exclusions"][stem]["Reason"],
            )
        elif item is None:
            status, reason = (
                "MISSING",
                "No FIA nuclei result for this source image",
            )
        else:
            mask, image = item
            row.update({key: image[key] for key in COUNTS})
            row["FIA_Mask_Name"] = mask
            row[LOW] = image[COUNT] < minimum
            try:
                row["FIA_Identity_Check"] = _identity(
                    row, image, bundle["calibration"].get(mask, {})
                )
                row["FN_Area_um2"] = _area_um2(row)
            except ValidationError as error:
                status, reason = "IDENTITY_OR_AREA_MISMATCH", str(error)
        row["FIA_Match_Status"] = status
        if status == "MATCHED":
            if row[COUNT] == 0:
                row["FN_Per_Nucleus_Reason"] = (
                    "Zero non-border nuclei; ratio undefined"
                )
            elif row["FN_Area_um2"] is None:
                row["FN_Per_Nucleus_Reason"] = (
                    "FN physical area/calibration unavailable"
                )
            else:
                row[RATIO] = row["FN_Area_um2"] / row[COUNT]
            reasons = [row[FN_REASON_COLUMN]] if row[FN_REASON_COLUMN] else []
            if row[LOW]:
                reasons.append(f"{COUNT} < {minimum}")
            row[FN_REASON_COLUMN] = "; ".join(reasons) or None
            row[FN_INCLUDED_FLAG] = not row[FN_LOW_FLAG] and not row[LOW]
            accepted.append(row)
            if row["FN_Per_Nucleus_Reason"]:
                log.event(
                    "WARNING",
                    "FN per nucleus unavailable",
                    f"{row['Image_ID']}: {row['FN_Per_Nucleus_Reason']}",
                )
            if not row[FN_INCLUDED_FLAG]:
                log.event(
                    "INFO",
                    "Image filter exclusion",
                    f"{row['Image_ID']}; well={row['Well']}; "
                    f"group={row['Group']}; FN%={row['FN_Area_Percent']:g}; "
                    f"nuclei={row[COUNT]}; "
                    f"FN cutoff={data['fn_threshold']:g}; "
                    f"nuclei cutoff={minimum}; {row[FN_REASON_COLUMN]}",
                    console=False,
                )
        else:
            row[FN_INCLUDED_FLAG] = False
            row[FN_REASON_COLUMN] = "; ".join(
                filter(None, [row[FN_REASON_COLUMN], reason])
            )
            log.event(
                "WARNING",
                "FIA image unavailable",
                f"{row['Image_ID']}: {reason}",
            )
        audit.append(row)
        matches.append(
            dict(
                zip(
                    MATCH_COLUMNS,
                    [row["Image_ID"], row["FIA_Mask_Name"], status, reason],
                )
            )
        )
    for stem, (mask, _) in lookup.items():
        if stem not in seen:
            matches.append(
                dict(
                    zip(
                        MATCH_COLUMNS,
                        [
                            stem,
                            mask,
                            "FIA_ONLY",
                            "Not present in the selected UMA collection",
                        ],
                    )
                )
            )
            log.event("INFO", "FIA extra image", mask, console=False)
    data["columns"].extend(ADDED_COLUMNS)
    data.update(
        fia_enabled=True,
        min_nuclei=minimum,
        fia_audit=audit,
        fia_matches=matches,
        fia_partial=len(accepted) != len(audit),
        rows=accepted,
        fia_run=str(bundle["run"]),
        fia_particle_size=bundle["particle_size"],
    )
    save_csv(directory / "fia_image_check.csv", MATCH_COLUMNS, matches)
    save_csv(directory / "fia_image_audit.csv", data["columns"], audit)
    if not accepted:
        raise ValidationError(
            "No usable matching FIA images; see fia_image_check.csv"
        )
    retained, excluded, groups, wells = filter_counts(
        accepted, data["group_order"], data["group_wells"]
    )
    data.update(
        retained_rows=retained,
        excluded_rows=excluded,
        group_filter_counts=groups,
        well_filter_counts=wells,
    )
    # Replace FN-only derived QC counts with the combined population.
    data["qc"] = [
        q
        for q in data["qc"]
        if q["Check"]
        not in {
            "Images in filtered data",
            "Excluded from filtered plots",
            "Empty filtered groups",
            "Single-image filtered groups",
        }
    ]
    for check in data["qc"]:
        if check["Check"] == "Generated plots":
            check.update(
                Value=8,
                Details="Four full-data and four combined-filter plots",
            )
        elif check["Check"] == "Images in full data":
            check.update(
                Value=len(accepted),
                Details=(
                    "FIA-matched images; all UMA inputs preserved "
                    "in FIA Image Audit"
                ),
            )
        elif check["Check"] == "Removed from full data":
            check.update(
                Value=len(audit) - len(accepted),
                Details="Unavailable FIA rows remain in FIA Image Audit",
            )
        elif (
            check["Check"] == "Overall validation status"
            and data["fia_partial"]
        ):
            check.update(
                Value="PARTIAL",
                Details=(
                    "Some UMA images lack usable FIA results; "
                    "see FIA Image Check"
                ),
            )
    data["qc"].extend(
        [
            {
                "Check": "FIA integration",
                "Value": bundle["run"].name,
                "Details": (
                    f"{len(accepted)}/{len(audit)} UMA images matched; "
                    f"{COUNT}; minimum {minimum}"
                ),
            },
            {
                "Check": "Images in filtered data",
                "Value": len(retained),
                "Details": (
                    "FN and nuclei thresholds are inclusive; both required"
                ),
            },
            {
                "Check": "Excluded from filtered plots",
                "Value": len(excluded),
                "Details": (
                    "Image-level filters; originals preserved "
                    "in FIA Image Audit"
                ),
            },
        ]
    )
    save_json(
        directory / "fia_selection.json",
        {
            "run": str(bundle["run"]),
            "particle_size_px2": bundle["particle_size"],
            "count_field": COUNT,
            "min_nuclei": minimum,
            "matched_images": len(accepted),
            "unmatched_images": len(audit) - len(accepted),
            "identity_checks": {
                r["Image_ID"]: r["FIA_Identity_Check"] for r in audit
            },
        },
    )
    log.event(
        "INFO",
        "Combined filters",
        f"{len(accepted)} matched images; {len(retained)} retained; "
        f"{len(excluded)} below FN and/or nuclei thresholds; "
        f"{len(audit) - len(accepted)} unavailable",
    )
