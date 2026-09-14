"""Console entry points. Parse arguments before importing image analysis tools."""

import argparse
from importlib.metadata import version
import logging
import sys


def _shutdown_imagej_workers():
    """Close ImageJ1's shared workers at the end of the standalone command.

    These non-daemon threads are outside the SciJava context, so disposing
    ImageJ alone leaves Python waiting for them at interpreter shutdown.
    Reusable analysis functions deliberately do not close this shared pool.
    """
    from scyjava import jimport, jvm_started

    if not jvm_started():
        return
    pool = jimport("ij.util.ThreadUtil").threadPoolExecutor
    seconds = jimport("java.util.concurrent.TimeUnit").SECONDS
    pool.shutdown()
    if not pool.awaitTermination(5, seconds):
        pool.shutdownNow()
        if not pool.awaitTermination(5, seconds):
            raise RuntimeError("ImageJ worker pool did not terminate")


def alignment():
    """Run the existing alignment workflow from the active environment."""
    parser = argparse.ArgumentParser(description="Fibronectin alignment analysis")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file containing folder_paths")
    parser.add_argument("-a", "--angle_value", type=float, default=15,
                        help="Alignment angle in degrees (default: 15)")
    args = parser.parse_args()
    from .alignment_analysis import main_fibronectin_processing
    main_fibronectin_processing(args.input, args.angle_value)


def thickness():
    """Run the existing thickness workflow from the active environment."""
    parser = argparse.ArgumentParser(description="Fibronectin thickness analysis")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {version('uma-tools')}")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file containing folder_paths")
    args = parser.parse_args()
    print(f"UMA-tools {version('uma-tools')} — thickness analysis", flush=True)
    from .thickness_analysis import main
    try:
        main(args.input)
    finally:
        analysis_failed = sys.exc_info()[0] is not None
        try:
            _shutdown_imagej_workers()
        except Exception:
            if not analysis_failed:
                raise
            # Keep the original analysis exception and its nonzero exit status.
            logging.exception("Could not close ImageJ workers after analysis failed.")


def area():
    """Run the standalone area command and preserve its process exit status."""
    from .area_analysis import main
    return main()


def collect_results():
    """Collect existing CSV results without importing image-analysis tools."""
    from .collect_results import main
    return main()


def report():
    """Create plots and Excel reports without starting image-analysis runtimes."""
    from .report import main
    return main()
