"""Command tests, with opt-in Java process-exit checks and no image data."""

import os
import subprocess
import sys
import textwrap
import types
import unittest
from unittest.mock import Mock, patch

from uma_tools import cli


class CommandTests(unittest.TestCase):
    def test_alignment_arguments(self):
        module = types.ModuleType("uma_tools.alignment_analysis")
        module.main_fibronectin_processing = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(
                sys, "argv", ["uma_alignment", "-i", "a b.json", "-a", "10"]
            ):
                cli.alignment()
        module.main_fibronectin_processing.assert_called_once_with(
            "a b.json", 10.0
        )

    def test_alignment_default_angle(self):
        module = types.ModuleType("uma_tools.alignment_analysis")
        module.main_fibronectin_processing = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(
                sys, "argv", ["uma_alignment", "-i", "input.json"]
            ):
                cli.alignment()
        module.main_fibronectin_processing.assert_called_once_with(
            "input.json", 15
        )

    def test_thickness_arguments(self):
        module = types.ModuleType("uma_tools.thickness_analysis")
        module.main = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(
                sys, "argv", ["uma_thickness", "-i", "a b.json"]
            ):
                cli.thickness()
        module.main.assert_called_once_with("a b.json")

    def test_area_preserves_analysis_exit_status(self):
        module = types.ModuleType("uma_tools.area_analysis")
        module.main = Mock(return_value=1)
        with patch.dict(sys.modules, {module.__name__: module}):
            self.assertEqual(cli.area(), 1)
        module.main.assert_called_once_with()

    def test_thickness_cleanup_preserves_analysis_error(self):
        module = types.ModuleType("uma_tools.thickness_analysis")
        module.main = Mock(side_effect=ValueError("original analysis error"))
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(
                sys, "argv", ["uma_thickness", "-i", "input.json"]
            ):
                with patch.object(
                    cli,
                    "_shutdown_imagej_workers",
                    side_effect=RuntimeError("cleanup error"),
                ) as cleanup:
                    with self.assertLogs(level="ERROR"):
                        with self.assertRaisesRegex(
                            ValueError, "original analysis error"
                        ):
                            cli.thickness()
        cleanup.assert_called_once_with()

    @unittest.skipUnless(
        os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
        "Set UMA_RUN_IMAGEJ_TESTS=1 to run Java process-exit checks",
    )
    def test_thickness_process_exits_after_java_work_on_success_and_failure(
        self,
    ):
        self.check_imagej_command_exit("thickness")

    @unittest.skipUnless(
        os.environ.get("UMA_RUN_IMAGEJ_TESTS") == "1",
        "Set UMA_RUN_IMAGEJ_TESTS=1 to run Java process-exit checks",
    )
    def test_alignment_process_exits_after_java_work_on_success_and_failure(
        self,
    ):
        self.check_imagej_command_exit("alignment")

    def check_imagej_command_exit(self, command):
        # A separate interpreter must exit naturally; shutting down workers in
        # the test runner would conceal the production command's exit bug.
        script = textwrap.dedent("""
            import os
            import sys
            import types
            import jpype
            from uma_tools import cli

            command = os.environ["UMA_TEST_COMMAND"]
            module = types.ModuleType(f"uma_tools.{command}_analysis")
            def main(*_):
                jar = os.environ.get("UMA_TEST_IMAGEJ_JAR")
                if jar:
                    jpype.startJVM("-Djava.awt.headless=true", classpath=[jar])
                else:
                    import imagej
                    ij = imagej.init("sc.fiji:fiji:2.14.0", mode="headless")
                pool = jpype.JClass("ij.util.ThreadUtil").threadPoolExecutor
                task = jpype.JProxy(
                    "java.lang.Runnable", dict(run=lambda: None)
                )
                pool.submit(task).get()
                assert pool.getPoolSize() > 0
                if not jar:
                    ij.dispose()
                print("ImageJ work finished", flush=True)
                if os.environ["UMA_TEST_ANALYSIS_FAILS"] == "1":
                    raise ValueError("deliberate analysis failure")

            entry_point = (
                "main_fibronectin_processing"
                if command == "alignment" else "main"
            )
            setattr(module, entry_point, main)
            sys.modules[module.__name__] = module
            sys.argv = [f"uma_{command}", "-i", "input.json"]
            getattr(cli, command)()
        """)
        for fails in (False, True):
            with self.subTest(analysis_fails=fails):
                env = dict(
                    os.environ,
                    UMA_TEST_ANALYSIS_FAILS=str(int(fails)),
                    UMA_TEST_COMMAND=command,
                )
                result = subprocess.run(
                    [sys.executable, "-c", script],
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                self.assertIn("ImageJ work finished", result.stdout)
                self.assertEqual(
                    result.returncode, 1 if fails else 0, result.stderr
                )
                if fails:
                    self.assertIn(
                        "ValueError: deliberate analysis failure",
                        result.stderr,
                    )

    def test_help_does_not_import_imagej(self):
        for function in ("alignment", "thickness", "area", "collect_results"):
            script = (
                "import sys; from uma_tools import cli; "
                "sys.argv=['assay','--help']; "
                "\ntry: getattr(cli, " + repr(function) + ")()"
                "\nexcept SystemExit as error:"
                "\n assert error.code == 0"
                "\nassert 'imagej' not in sys.modules"
                "\nassert 'matplotlib.pyplot' not in sys.modules"
            )
            result = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_thickness_version_does_not_start_imagej(self):
        script = (
            "import sys; from uma_tools import cli; "
            "sys.argv=['uma_thickness','--version']; "
            "\ntry: cli.thickness()"
            "\nexcept SystemExit as error:"
            "\n assert error.code == 0"
            "\nassert 'imagej' not in sys.modules"
        )
        result = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("uma_thickness", result.stdout)

    def test_installed_commands_outside_repository(self):
        for command in (
            "uma_alignment",
            "uma_thickness",
            "area_analysis",
            "uma_collect_results",
        ):
            executable = os.path.join(os.path.dirname(sys.executable), command)
            result = subprocess.run(
                [executable, "--help"],
                cwd="/tmp",
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--input", result.stdout)
            result = subprocess.run(
                [executable],
                cwd="/tmp",
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(result.returncode, 2, result.stderr)


if __name__ == "__main__":
    unittest.main()
