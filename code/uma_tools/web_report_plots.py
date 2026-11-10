"""Interactive geometry from saved UMA plots; no scientific recalculation."""

from .web_report_data import numeric, rows


def charts(run):
    """Use saved boxes, points, colors, brackets, limits and populations."""
    if run["kind"] != "Report" or run["status"] not in ("SUCCESS", "PARTIAL"):
        return []
    if not run["integrity_ok"]:
        return []
    manifest = run["documents"].get("plot_manifest.json", [])
    if not manifest:
        run["warnings"].append(
            "No saved plot manifest; tables remain available"
        )
        return []
    result = []
    all_rows = rows(run, "plot_data")
    for saved in manifest:
        plot = saved["plot_id"]
        if saved.get("name") not in (
            "Alignment",
            "Area",
            "Fibronectin",
            "FN per nucleus",
        ):
            continue
        points = [r for r in all_rows if r.get("Plot") == plot]
        ids = [r["Image_ID"] for r in points]
        if (
            len(ids) != len(set(ids))
            or set(ids) != set(saved["plotted_image_ids"])
            or len(ids) != saved["point_count"]
        ):
            run["warnings"].append(f"{plot}: saved plot population mismatch")
            continue
        if any(
            numeric(r.get("Value")) is None
            or numeric(r.get("Panel_X")) is None
            for r in points
        ):
            run["warnings"].append(f"{plot}: invalid saved point coordinates")
            continue
        result.append({**saved, "points": points})
    return result
