"""Run unchanged imaging engines on scheme-selected marker/channel subsets."""

from pathlib import Path

from .config import read_config
from .fiber_routes import preflight, routes, use_route
from .progress import (
    CANCELLATIONS,
    console,
    folder_error,
    folder_scope,
    outcome,
    register_sources,
)


def run_images(step, args):
    """Preflight the complete batch, then share one ImageJ context."""
    if step == "area":
        from . import area_analysis as area
        from .area_imagej import DEFAULT_THRESHOLD_LOWER, threshold_settings

        folders, input_json = area.read_source_folders(args)
        bounds = args.threshold or [DEFAULT_THRESHOLD_LOWER]
        if not 1 <= len(bounds) <= 2:
            raise ValueError("Use -t LOWER or -t LOWER UPPER")
        limits = (
            None
            if args.projections_only
            else threshold_settings(
                bounds[0], bounds[1] if len(bounds) == 2 else None
            )
        )
    else:
        input_json = Path(args.input).expanduser().resolve()
        folders = read_config(input_json)
    register_sources(folders, input_json)
    prepared = preflight(folders, args)
    holder, outputs, ij = [None], [], None
    failed = succeeded = 0
    try:
        for item, design in prepared:
            with folder_scope(item.source):
                if item.error:
                    folder_error(item.source, item.error)
                    failed += 1
                    continue
                plans = routes(
                    item.source, design, getattr(args, "marker", None)
                )
                if not plans:
                    outcome(
                        "NO_INPUT",
                        "No assigned images for requested Fiber marker",
                    )
                    failed += 1
                    continue
                folder_failed = False
                folder_usable = False
                states = []
                for route in plans:
                    try:
                        with use_route(route):
                            if step == "area":
                                result = area.process_folder(
                                    item.source,
                                    route.channel,
                                    args.projection,
                                    limits,
                                    holder,
                                    input_json,
                                    outputs,
                                ).value
                            else:
                                from .imagej import initialize_imagej

                                if ij is None and any(
                                    r["Status"] == "SELECTED"
                                    for r in route.inventory()
                                ):
                                    ij = initialize_imagej()
                                if step == "alignment":
                                    from .alignment_analysis import (
                                        process_folder,
                                    )

                                    result = process_folder(
                                        str(item.source),
                                        route.channel,
                                        args.angle_value,
                                        500,
                                        500,
                                        ij,
                                    )
                                else:
                                    from .thickness_analysis import (
                                        import_java_classes,
                                        process_single_folder,
                                    )

                                    classes = (
                                        import_java_classes()
                                        if ij is not None
                                        else (None,) * 6
                                    )
                                    IJ, _, windows, table, duplicator, _ = (
                                        classes
                                    )
                                    result = process_single_folder(
                                        ij,
                                        IJ,
                                        windows,
                                        duplicator,
                                        table,
                                        str(item.source),
                                        None,
                                        route.channel,
                                    )
                            failed += result != "SUCCESS"
                            states.append(result)
                            succeeded += result == "SUCCESS"
                            folder_failed |= result != "SUCCESS"
                            folder_usable |= result in ("SUCCESS", "PARTIAL")
                    except CANCELLATIONS:
                        raise
                    except Exception as error:
                        folder_error(item.source, error)
                        failed += 1
                        folder_failed = True
                        states.append("ERROR")
                if folder_failed:
                    outcome(
                        "NO_INPUT"
                        if states and set(states) == {"NO_INPUT"}
                        else "PARTIAL"
                        if folder_usable
                        else "FAILED",
                        f"{item.source}: some Fiber routes failed or are "
                        f"incomplete",
                    )
    finally:
        if step == "area":
            if not area.finish_imagej(holder, outputs):
                failed += 1
        elif ij is not None:
            ij.dispose()
    console(
        f"{step}: {succeeded} successful marker/channel runs; {failed} "
        f"incomplete or failed."
    )
    return int(bool(failed))
