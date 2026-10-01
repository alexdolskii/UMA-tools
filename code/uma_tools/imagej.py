"""Initialize Fiji and close its workers at command boundaries."""

from contextlib import contextmanager
from typing import Any

from .progress import current_folder, current_session, phase

FIJI_ENDPOINT = "sc.fiji:fiji:2.14.0"


class ImageJInitializationError(Exception):
    """Fiji could not initialize its headless context."""


def initialize_imagej() -> Any:
    """Return the established headless Fiji context."""
    phase("Initializing ImageJ")
    try:
        import imagej

        context = imagej.init(FIJI_ENDPOINT, mode="headless")
    except Exception as error:
        raise ImageJInitializationError(
            f"Failed to initialize ImageJ: {error}"
        ) from error
    phase("ImageJ initialization completed")
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


class _BioformatsAppender:
    """Logback adapter; keep Java messages out of the progress line."""

    def __init__(self, session, folder, jclass):
        self.session, self.folder, self.jclass = session, folder, jclass
        self.name = "UMAProgress"
        self.started = True
        self.context = None
        self.errors = []

    def doAppend(self, event):
        severity = event.getLevel().toInt()
        level = (
            "ERROR"
            if severity >= 40000
            else "WARNING"
            if severity >= 30000
            else "INFO"
        )
        message = str(event.getFormattedMessage())
        throwable = event.getThrowableProxy()
        if throwable is not None:
            utility = self.jclass(
                "ch.qos.logback.classic.spi.ThrowableProxyUtil"
            )
            message += "\n" + str(utility.asString(throwable))
        try:
            self.session.event(
                level, "Bio-Formats", message, folder=self.folder
            )
            if level != "INFO":
                self.session.progress.message(
                    f"{level}: {message.splitlines()[0]}"
                )
        except OSError as error:
            # Java logging swallows appender exceptions. Surface them in
            # Python at the end of the image operation instead.
            self.errors.append(error)

    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    def start(self):
        self.started = True

    def stop(self):
        self.started = False

    def isStarted(self):
        return self.started

    def getContext(self):
        return self.context

    def setContext(self, context):
        self.context = context

    def addStatus(self, *_):
        pass

    addInfo = addWarn = addError = addFilter = clearAllFilters = addStatus

    def getCopyOfAttachedFiltersList(self):
        return self.jclass("java.util.ArrayList")()

    def getFilterChainDecision(self, _):
        return self.jclass("ch.qos.logback.core.spi.FilterReply").NEUTRAL


@contextmanager
def bioformats_log():
    """
    Restore the original Java logger after an image, even on failure.
    """
    session = current_session()
    if session is None:
        yield
        return
    import jpype

    if not jpype.isJVMStarted():
        yield
        return
    logger = proxy = adapter = None
    configured = False
    previous_appenders = []
    try:
        logger = jpype.JClass("org.slf4j.LoggerFactory").getLogger(
            "loci.formats"
        )
        if not isinstance(
            logger, jpype.JClass("ch.qos.logback.classic.Logger")
        ):
            raise TypeError("Bio-Formats does not use Logback")
        previous_level = logger.getLevel()
        previous_additive = logger.isAdditive()
        previous_appenders = list(logger.iteratorForAppenders())
        adapter = _BioformatsAppender(session, current_folder(), jpype.JClass)
        proxy = jpype.JProxy("ch.qos.logback.core.Appender", inst=adapter)
        configured = True
        for appender in previous_appenders:
            logger.detachAppender(appender)
        logger.addAppender(proxy)
        logger.setAdditive(False)
        logger.setLevel(jpype.JClass("ch.qos.logback.classic.Level").INFO)
    except Exception as error:
        if configured:
            logger.detachAppender(proxy)
            logger.setLevel(previous_level)
            logger.setAdditive(previous_additive)
            for appender in previous_appenders:
                logger.addAppender(appender)
            configured = False
        if not getattr(session, "bioformats_warning", False):
            session.event(
                "WARNING",
                "Bio-Formats",
                "Compact logging unavailable; normal output retained: "
                f"{error}",
                announce=True,
            )
            session.bioformats_warning = True
    try:
        yield
        if adapter is not None and adapter.errors:
            raise OSError(f"Bio-Formats log write failed: {adapter.errors[0]}")
    finally:
        if configured:
            logger.detachAppender(proxy)
            logger.setLevel(previous_level)
            logger.setAdditive(previous_additive)
            for appender in previous_appenders:
                logger.addAppender(appender)
