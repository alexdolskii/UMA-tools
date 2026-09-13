"""Tests for the installed commands; no JVM or image data is needed."""

import os
import subprocess
import sys
import types
import unittest
from unittest.mock import Mock, patch

from uma_tools import cli


class CommandTests(unittest.TestCase):
    def test_alignment_arguments(self):
        module = types.ModuleType("uma_tools.alignment_analysis")
        module.main_fibronectin_processing = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(sys, "argv", ["uma_alignment", "-i", "a b.json", "-a", "10"]):
                cli.alignment()
        module.main_fibronectin_processing.assert_called_once_with("a b.json", 10.0)

    def test_alignment_default_angle(self):
        module = types.ModuleType("uma_tools.alignment_analysis")
        module.main_fibronectin_processing = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(sys, "argv", ["uma_alignment", "-i", "input.json"]):
                cli.alignment()
        module.main_fibronectin_processing.assert_called_once_with("input.json", 15)

    def test_thickness_arguments(self):
        module = types.ModuleType("uma_tools.thickness_analysis")
        module.main = Mock()
        with patch.dict(sys.modules, {module.__name__: module}):
            with patch.object(sys, "argv", ["uma_thickness", "-i", "a b.json"]):
                cli.thickness()
        module.main.assert_called_once_with("a b.json")

    def test_help_does_not_import_imagej(self):
        for function in ("alignment", "thickness"):
            script = (
                "import sys; from uma_tools import cli; sys.argv=['assay','--help']; "
                "\ntry: getattr(cli, " + repr(function) + ")()"
                "\nexcept SystemExit as error:"
                "\n assert error.code == 0"
                "\nassert 'imagej' not in sys.modules"
                "\nassert 'matplotlib.pyplot' not in sys.modules"
            )
            result = subprocess.run([sys.executable, "-c", script],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_installed_commands_outside_repository(self):
        for command in ("uma_alignment", "uma_thickness"):
            executable = os.path.join(os.path.dirname(sys.executable), command)
            result = subprocess.run([executable, "--help"], cwd="/tmp",
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--input", result.stdout)
            result = subprocess.run([executable], cwd="/tmp",
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 2, result.stderr)


if __name__ == "__main__":
    unittest.main()
