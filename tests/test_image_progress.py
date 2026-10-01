"""Stage progress must reflect work without finalizing measurements early."""

import io
import json
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import Mock, patch

from uma_tools import thickness_analysis as thickness
from uma_tools.image_run import ImageRun
from uma_tools.progress import (
    CANCELLATIONS,
    CommandSession,
    CompactProgress,
    folder_scope,
    phase,
)
from uma_tools.run import RunLog


class TerminalBuffer(io.StringIO):
    def isatty(self):
        return True


class ImageProgressTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        self.sources = [root / "plate1", root / "plate2"]
        for source in self.sources:
            source.mkdir()
        self.terminal = TerminalBuffer()
        self.session = CommandSession("alignment", self.sources)
        self.session.progress = CompactProgress(self.terminal)
        contexts = ExitStack()
        self.addCleanup(contexts.close)
        contexts.enter_context(
            patch(
                "uma_tools.progress.shutil.get_terminal_size",
                return_value=Mock(columns=240),
            )
        )
        contexts.enter_context(self.session)

    def make_run(self, source, names):
        output = source / "uma_assay" / "test_results"
        output.mkdir()
        return ImageRun(output, source, names, {})

    def assert_progress(
        self, source, stage, done, total, name, failed=0, *, label=None
    ):
        expected = f"{source.name}: {stage}: {done}/{total} finished"
        if failed:
            expected += f" | {failed} failed"
        expected += f" | {name}"
        self.assertEqual(
            self.session.progress.label if label is None else label, expected
        )
        self.assertIn(expected, self.terminal.getvalue())

    def test_alignment_updates_each_stage_after_every_image(self):
        source = self.sources[0]
        names = ("one.nd2", "two.nd2")
        run = self.make_run(source, names)
        with folder_scope(source), run:
            for stage in ("Projection", "Orientation", "Summary"):
                for index, name in enumerate(names):
                    with run.attempt(name, stage, final=stage == "Summary"):
                        during = self.session.progress.label
                    self.assert_progress(
                        source, stage, index, 2, name, label=during
                    )
                    self.assert_progress(source, stage, index + 1, 2, name)
                    state = json.loads(
                        (run.output / "run_status.json").read_text()
                    )
                    self.assertEqual(
                        state["processed_images"],
                        index + 1 if stage == "Summary" else 0,
                    )
            self.assertEqual(run.finish(), "SUCCESS")
        journal = (source / "uma_assay/UMA_Logs/1_alignment.log").read_text()
        for stage in ("Projection", "Orientation", "Summary"):
            self.assertIn(f"{stage}: 2/2 finished", journal)

    def test_failed_images_finish_the_stage_and_leave_later_denominators(self):
        source = self.sources[0]
        names = (
            "same.nd2",
            "same.tif",
            "p.nd2",
            "o.nd2",
            "s.nd2",
            "good.nd2",
        )
        run = self.make_run(source, names)
        with folder_scope(source), run:
            run.reject_collisions()
            for stage, bad, total in (
                ("Projection", "p.nd2", 4),
                ("Orientation", "o.nd2", 3),
                ("Summary", "s.nd2", 2),
            ):
                selected = [name for name in names if run.eligible(name)]
                self.assertEqual(len(selected), total)
                for index, name in enumerate(selected):
                    with run.attempt(name, stage, final=stage == "Summary"):
                        during = self.session.progress.label
                        if name == bad:
                            raise ValueError("Deliberate unreadable image")
                    self.assert_progress(
                        source,
                        stage,
                        index,
                        total,
                        name,
                        int(index > 0),
                        label=during,
                    )
                    self.assert_progress(
                        source, stage, index + 1, total, name, failed=1
                    )
            self.assertEqual(run.finish(), "PARTIAL")
        self.assertEqual(run.status["processed_images"], 1)
        self.assertEqual(run.status["failed_images"], 5)
        self.assertEqual(run.status["unprocessed_images"], 0)

    def test_cancellation_never_advances_the_counter(self):
        source = self.sources[0]
        for index, cancellation in enumerate(CANCELLATIONS):
            with self.subTest(cancellation=cancellation.__name__):
                output = source / "uma_assay" / f"cancel_{index}"
                output.mkdir()
                run = ImageRun(output, source, ["one.nd2", "two.nd2"], {})
                with self.assertRaises(cancellation):
                    with folder_scope(source), run:
                        with run.attempt("one.nd2", "Projection"):
                            pass
                        with run.attempt("two.nd2", "Projection"):
                            raise cancellation()
                self.assert_progress(source, "Projection", 1, 2, "two.nd2")
                self.assertEqual(run.status["status"], "CANCELLED")
                self.assertEqual(run.status["processed_images"], 0)
                self.assertEqual(run.status["failed_images"], 0)

    def test_single_stage_counters_reset_in_the_next_folder(self):
        for stage in ("Thickness", "Projection and area"):
            with self.subTest(stage=stage):
                for source, total in zip(self.sources, (2, 1)):
                    output = source / "uma_assay" / stage
                    output.mkdir()
                    names = [f"image{i}.tif" for i in range(total)]
                    run = ImageRun(output, source, names, {})
                    with folder_scope(source), run:
                        for index, name in enumerate(names):
                            with run.attempt(name, stage, final=True):
                                during = self.session.progress.label
                            self.assert_progress(
                                source,
                                stage,
                                index,
                                total,
                                name,
                                label=during,
                            )
                            self.assert_progress(
                                source, stage, index + 1, total, name
                            )
                        self.assertEqual(run.finish(), "SUCCESS")

    def test_operations_keep_the_counter_and_current_filename(self):
        source = self.sources[0]
        names = ("one.nd2", "two.nd2")
        run = self.make_run(source, names)
        observed = []
        with folder_scope(source), run:
            for index, name in enumerate(names):
                with run.attempt(name, "Thickness", final=True):
                    for operation in (
                        "Opening image...",
                        "Running Local Thickness...",
                        "Closed all images.\n",
                    ):
                        phase(operation)
                        observed.append(
                            (
                                index,
                                name,
                                operation.strip(),
                                self.session.progress.label,
                            )
                        )
            run.finish()
            phase("Closing ImageJ context")
            self.assertEqual(
                self.session.progress.label, "plate1: Closing ImageJ context"
            )
        for index, name, operation, label in observed:
            expected = (
                f"plate1: Thickness: {index}/2 finished | {operation} | {name}"
            )
            self.assertEqual(label, expected)
            self.assertIn(expected, self.terminal.getvalue())
        journal = (source / "uma_assay/UMA_Logs/1_alignment.log").read_text()
        self.assertIn(
            "Thickness: 1/2 finished | Opening image... | two.nd2", journal
        )

    def test_real_thickness_filter_messages_do_not_replace_progress(self):
        source = self.sources[0]
        run = self.make_run(source, ["one.nd2"])
        calls = []
        ij = Mock()
        ij.run.side_effect = lambda *args: calls.append(
            (args, self.session.progress.label)
        )
        image = object()
        with folder_scope(source), run:
            with run.attempt("one.nd2", "Thickness", final=True):
                thickness.create_thickness_mask(ij, image)
            run.finish()
        expected_operations = (
            "Applying Maximum filter...",
            "Applying Gaussian Blur...",
            "Subtracting background...",
            "Applying threshold...",
            "Applying threshold...",
        )
        self.assertEqual(len(calls), len(expected_operations))
        for (_, label), operation in zip(calls, expected_operations):
            self.assertEqual(
                label,
                f"plate1: Thickness: 0/1 finished | {operation} | one.nd2",
            )

    def test_failure_and_cancellation_clear_the_image_activity_context(self):
        source = self.sources[0]
        run = self.make_run(source, ["bad.nd2", "next.nd2"])
        with folder_scope(source), run:
            with run.attempt("bad.nd2", "Thickness", final=True):
                phase("Opening image...")
                raise ValueError("Deliberate unreadable image")
            phase("Between images")
            self.assertEqual(
                self.session.progress.label, "plate1: Between images"
            )
            with run.attempt("next.nd2", "Thickness", final=True):
                phase("Opening image...")
                observed = self.session.progress.label
            self.assertEqual(
                observed,
                "plate1: Thickness: 1/2 finished | 1 failed | "
                "Opening image... | next.nd2",
            )
            self.assertEqual(run.finish(), "PARTIAL")
        for index, cancellation in enumerate(CANCELLATIONS):
            output = source / "uma_assay" / f"cancel_activity_{index}"
            output.mkdir()
            with self.assertRaises(cancellation):
                with (
                    folder_scope(source),
                    ImageRun(output, source, ["cancel.nd2"], {}) as run,
                ):
                    with run.attempt("cancel.nd2", "Thickness"):
                        phase("Opening image...")
                        raise cancellation()
            with folder_scope(self.sources[1]):
                phase("Initializing ImageJ")
                self.assertEqual(
                    self.session.progress.label, "plate2: Initializing ImageJ"
                )

    def test_area_events_also_keep_the_counter_without_duplicate_log_entries(
        self,
    ):
        source = self.sources[0]
        run = self.make_run(source, ["one.nd2"])
        with folder_scope(source), run:
            log = RunLog(run.output)
            try:
                with run.attempt("one.nd2", "Projection and area", final=True):
                    log.event("PASS", "Projection", "SUM32 verified")
                    observed = self.session.progress.label
                    log.event("WARNING", "Area", "Example warning")
                    warning = self.session.progress.label
                run.finish()
            finally:
                log.close()
        self.assertEqual(
            observed,
            "plate1: Projection and area: 0/1 finished | "
            "Projection: SUM32 verified | one.nd2",
        )
        self.assertEqual(
            warning,
            "plate1: Projection and area: 0/1 finished | "
            "Area: Example warning | one.nd2",
        )
        journal = (source / "uma_assay/UMA_Logs/1_alignment.log").read_text()
        self.assertEqual(journal.count("SUM32 verified"), 1)
