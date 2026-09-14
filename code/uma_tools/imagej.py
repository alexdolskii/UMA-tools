"""Initialize Fiji and close its workers at command boundaries."""

from typing import Any

FIJI_ENDPOINT = "sc.fiji:fiji:2.14.0"


class ImageJInitializationError(Exception):
    """Fiji could not initialize its headless context."""


def initialize_imagej() -> Any:
    """Return the established headless Fiji context."""
    print("Initializing ImageJ...", flush=True)
    try:
        import imagej

        context = imagej.init(FIJI_ENDPOINT, mode="headless")
    except Exception as error:
        raise ImageJInitializationError(
            f"Failed to initialize ImageJ: {error}"
        ) from error
    print("ImageJ initialization completed.", flush=True)
    return context


def shutdown_imagej_workers() -> None:
    """Close ImageJ1 workers when a standalone command finishes.

    Context disposal does not stop these non-daemon threads. Close them
    at the command boundary, never between source folders or inside
    reusable image-analysis functions.
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
