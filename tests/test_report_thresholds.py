"""Distinguish saved mask thresholds from report coverage filtering."""

from pathlib import Path
from unittest.mock import patch

from test_plot_style import render, two_panels
from test_report import ReportFixture

from uma_tools.report import parse_args
from uma_tools.report_schema import ValidationError
from uma_tools.report_tables import prepare_display_tables

FLOAT32_MAX = float.fromhex("0x1.fffffep+127")


class MaskThresholdReportTests(ReportFixture):
    @staticmethod
    def thresholds(**overrides):
        values = {
            "Threshold_Lower": 2000,
            "Threshold_Upper": "",
            "Effective_Threshold_Upper": FLOAT32_MAX,
            "Threshold_Units": "Raw projection intensity",
        }
        values.update(overrides)
        return values

    def test_mask_intensity_is_not_the_coverage_cutoff(self):
        paths = self.inputs(
            percentages=[19, 20], mask_thresholds=self.thresholds()
        )
        originals = {key: path.read_bytes() for key, path in paths.items()}
        data = self.merge(paths)
        self.assertEqual(parse_args(["-i", "unused.json"]).fn_threshold, 20)
        self.assertEqual(data["fn_threshold"], 20)
        self.assertEqual(len(data["retained_rows"]), 1)
        self.assertEqual(len(data["excluded_rows"]), 1)
        settings = data["fn_mask_settings"]
        self.assertEqual(
            settings["caption"],
            "FN mask: SUM projection; raw intensity "
            "[2000, float32 max] (inclusive).",
        )
        self.assertTrue(settings["complete"])
        self.assertEqual(settings["profiles"][0]["images"], 2)
        self.assertEqual(settings["profiles"][0]["upper"], FLOAT32_MAX)
        for key, path in paths.items():
            self.assertEqual(path.read_bytes(), originals[key])

    def test_recorded_bounds_and_projection_override_no_report_defaults(self):
        paths = self.inputs(
            mask_thresholds=self.thresholds(
                Threshold_Lower=1234.5,
                Threshold_Upper=45000,
                Effective_Threshold_Upper=45000,
                Requested_Threshold_Lower=50,
            )
        )
        for index in range(2):
            self.mutate_csv(
                paths["fibronectin"], "Projection_Method", "MAX", index
            )
        data = self.merge(paths)
        self.assertEqual(
            data["fn_mask_settings"]["caption"],
            "FN mask: MAX projection; raw intensity "
            "[1234.5, 45000] (inclusive).",
        )

    def test_legacy_and_partial_metadata_do_not_invent_missing_thresholds(
        self,
    ):
        paths = self.inputs()
        data = self.merge(paths)
        self.assertIn(
            "intensity thresholds not recorded",
            data["fn_mask_settings"]["caption"],
        )
        self.assertFalse(data["fn_mask_settings"]["complete"])
        partial = self.inputs(
            combined=self.make_combined("20260914_120001_000001"),
            mask_thresholds={"Threshold_Lower": 3000},
        )
        data = self.merge(partial)
        caption = data["fn_mask_settings"]["caption"]
        self.assertIn("lower=3000, upper=not recorded", caption)
        self.assertIn("units not recorded", caption)
        self.assertNotIn("float32 max", caption)
        self.assertNotIn("2000", caption)
        self.assertEqual(len(data["rows"]), 2)
        self.assertTrue(any(e["Level"] == "WARNING" for e in self.log.events))

    def test_mixed_thresholds_are_reported_with_per_profile_image_counts(self):
        paths = self.inputs(mask_thresholds=self.thresholds())
        self.mutate_csv(paths["fibronectin"], "Threshold_Lower", 3000)
        data = self.merge(paths)
        settings = data["fn_mask_settings"]
        self.assertTrue(settings["mixed"])
        self.assertIn(
            "mixed intensity/projection settings (2)", settings["caption"]
        )
        self.assertEqual(
            {p["lower"]: p["images"] for p in settings["profiles"]},
            {2000: 1, 3000: 1},
        )
        prepare_display_tables(data, [], {}, [])
        tables = data["display_tables"]
        overview = {
            row["Item"]: row["Value"] for row in tables["Overview"]["rows"]
        }
        descriptions = " ".join(str(value) for value in overview.values())
        self.assertIn("raw intensity [2000, float32 max]", descriptions)
        self.assertIn("raw intensity [3000, float32 max]", descriptions)
        self.assertEqual(len(data["retained_rows"]), 2)

    def test_invalid_recorded_bounds_fail_instead_of_mislabelling_plots(self):
        paths = self.inputs(mask_thresholds=self.thresholds())
        for value in ("NaN", "inf", "not numeric", "-1", "4e38"):
            with self.subTest(value=value):
                self.mutate_csv(paths["fibronectin"], "Threshold_Lower", value)
                with self.assertRaisesRegex(
                    ValidationError, "FN mask.*source row"
                ):
                    self.merge(paths)

    def test_plot_exports_show_both_thresholds_without_changing_observations(
        self,
    ):
        # Load pyplot before replacing the Figure.savefig function.
        import matplotlib.pyplot  # noqa: F401

        data = two_panels()
        data["source_name"] = (
            "Synthetic plate with distinct mask and coverage thresholds"
        )
        data["run_id"] = "20261001_200000"
        paths = self.inputs(mask_thresholds=self.thresholds())
        mask_settings = self.merge(paths)["fn_mask_settings"]
        with patch("matplotlib.figure.Figure.savefig"):
            baselines = [
                render(data, self.root, mode) for mode in (False, True)
            ]
        data["fn_mask_settings"] = mask_settings

        from matplotlib.figure import Figure

        savefig = Figure.savefig
        captured = []

        def inspect(figure, *args, **kwargs):
            figure.canvas.draw()
            footer = next(
                ax for ax in figure.axes if ax.get_label() == "footer"
            )
            caption, source = footer.texts
            text = caption.get_text().replace("\n", " ")
            self.assertIn("raw intensity [2000, float32 max]", text)
            self.assertIn("Coverage filter: FN ≥ 20%", text)
            caption_box = caption.get_window_extent()
            source_box = source.get_window_extent()
            self.assertGreater(caption_box.y0, source_box.y1)
            for bounds in (caption_box, source_box):
                self.assertGreaterEqual(bounds.x0, figure.bbox.x0)
                self.assertGreaterEqual(bounds.y0, figure.bbox.y0)
                self.assertLessEqual(bounds.x1, figure.bbox.x1)
                self.assertLessEqual(bounds.y1, figure.bbox.y1)
            captured.append(text)
            return savefig(figure, *args, **kwargs)

        with patch.object(Figure, "savefig", inspect):
            for filtered, baseline in zip((False, True), baselines):
                result = render(data, self.root, filtered, "both")
                for key in (
                    "boxes",
                    "plotted_image_ids",
                    "red_outline_image_ids",
                    "statistical_comparisons",
                    "technical_replicate_colors",
                    "rendered_x_positions",
                    "y_min",
                    "y_max",
                ):
                    self.assertEqual(result[key], baseline[key], key)
                self.assertTrue(Path(result["png_file"]).is_file())
                self.assertTrue(
                    Path(result["pdf_file"]).read_bytes().startswith(b"%PDF")
                )
        self.assertEqual(len(captured), 4)
