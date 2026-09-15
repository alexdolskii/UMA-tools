"""Verify plot population, axes and supplied statistical annotations."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

from uma_tools import report_plots as plots
from uma_tools.report_schema import FN_LOW_FLAG, FN_METRIC


def report_data(statistics=None):
    """Small plate with retained, excluded and empty conditions."""
    records = [
        ("Control", "B02", "control_1", 50),
        ("Control", "B03", "control_2", 60),
        ("Treatment", "C02", "treatment_1", 70),
        ("Treatment", "C03", "treatment_2", 10),
        ("Empty", "D02", "empty_1", 5),
    ]
    rows = []
    for group, well, image_id, fn in records:
        rows.append(
            {
                "Group": group,
                "Well": well,
                "Image_ID": image_id,
                FN_METRIC: fn,
                FN_LOW_FLAG: fn < 20,
                "Alignment": fn / 2,
                "Area (µm²)": fn * 10,
                "StdDev (µm)": fn / 20,
                "Min (µm)": fn / 10,
                "Max (µm)": fn / 5,
                "Median (µm)": fn / 8,
            }
        )
    return {
        "fn_threshold": 20,
        "metric": "Alignment",
        "angle_label": "15",
        "group_order": ["Control", "Treatment", "Empty"],
        "group_wells": {
            "Control": ["B02", "B03"],
            "Treatment": ["C02", "C03"],
            "Empty": ["D02"],
        },
        "rows": rows,
        "retained_rows": [row for row in rows if not row[FN_LOW_FLAG]],
        "excluded_rows": [row for row in rows if row[FN_LOW_FLAG]],
        "statistics": statistics,
    }


def comparison(treatment="Treatment", **overrides):
    """Supply engine output; renderer must not reinterpret p-values."""
    result = {
        "Comparison_Block": "block_1",
        "Color_Code": "theme:6;tint:0",
        "Control": "Control",
        "Treatment": treatment,
        "Metric": FN_METRIC,
        "Stats_Unit": "well",
        "Status": "Tested",
        "Significance": "**",
        "P_Holm": 0.4,
    }
    result.update(overrides)
    return result


class StatisticalPlotsTests(unittest.TestCase):
    def setUp(self):
        import matplotlib

        matplotlib.use("Agg")
        # Load pyplot before replacing Figure.savefig's function object.
        import matplotlib.pyplot  # noqa: F401

        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.log = Mock()

    def render(self, data, filtered=True, spec_index=0):
        captured = []

        def capture(figure, *args, **kwargs):
            captured.append(figure)

        spec = plots._plot_specs(data)[spec_index]
        upper = {
            field: (
                100.0
                if name in ("Alignment", "Fibronectin")
                else max(row[field] for row in data["rows"]) * 1.12
            )
            for name, field, _, _ in plots._plot_specs(data)
        }
        with patch("matplotlib.figure.Figure.savefig", capture):
            manifest = plots._render_plot(
                data,
                self.directory,
                spec,
                filtered,
                plots.replicate_colors(2),
                upper,
                plots._point_positions(data, 2),
                self.log,
            )
        return manifest, captured[0]

    def test_disabled_report_includes_filtered_fn_without_annotations(self):
        with patch("matplotlib.figure.Figure.savefig"):
            manifest = plots.create_plots(
                report_data(), self.directory / "plots", self.log
            )
        self.assertEqual(len(manifest), 14)
        self.assertEqual(manifest[7]["sheet"], "Fibronectin Filtered")
        self.assertEqual(manifest[7]["point_count"], 3)
        self.assertEqual(manifest[7]["red_outline_count"], 0)
        for item in manifest:
            self.assertIsNone(item["statistics_unit"])
            self.assertEqual(item["statistical_comparisons"], [])
            self.assertEqual(
                item["statistics_note"], "Statistical tests disabled."
            )

    def test_statistics_preserve_filtered_points_colors_and_axes(self):
        baseline, old_figure = self.render(report_data())
        statistics = {
            "unit": "well",
            "comparisons": [
                comparison(),
                comparison("Empty", Status="Not tested", Significance=""),
            ],
        }
        manifest, figure = self.render(report_data(statistics))
        stable_fields = (
            "y_min",
            "y_max",
            "group_order",
            "group_counts",
            "group_well_counts",
            "plotted_image_ids",
            "red_outline_image_ids",
            "x_positions",
            "technical_replicate_colors",
            "boxes",
        )
        for field in stable_fields:
            self.assertEqual(baseline[field], manifest[field], field)
        self.assertEqual(len(figure.axes), 2)
        self.assertEqual(figure.axes[-1].get_ylim(), (0, 100))
        self.assertEqual(
            figure.axes[-1].get_xlim(), old_figure.axes[-1].get_xlim()
        )
        for old, new in zip(
            old_figure.axes[-1].collections, figure.axes[-1].collections
        ):
            self.assertEqual(
                old.get_facecolors().tolist(), new.get_facecolors().tolist()
            )
            self.assertEqual(
                old.get_offsets().tolist(), new.get_offsets().tolist()
            )
        self.assertEqual(
            [text.get_text() for text in figure.axes[0].texts],
            ["**", "Not tested"],
        )
        self.assertIn("one mean per well", manifest["statistics_note"])
        self.assertTrue(
            all(
                "n_images=" in label.get_text()
                and "n_wells=" in label.get_text()
                for label in figure.axes[-1].get_xticklabels()
            )
        )

    def test_full_view_never_gets_filtered_test_annotations(self):
        statistics = {"unit": "image", "comparisons": [comparison()]}
        manifest, figure = self.render(report_data(statistics), False)
        self.assertEqual(manifest["statistical_comparisons"], [])
        self.assertEqual(manifest["point_count"], 5)
        self.assertEqual(manifest["red_outline_count"], 2)
        self.assertEqual(len(figure.axes), 1)
        self.assertIn(
            "No tests on this full-data", manifest["statistics_note"]
        )

    def test_metric_selection_ns_and_image_dependence_note(self):
        statistics = {
            "unit": "image",
            "comparisons": [
                comparison(Significance="ns"),
                comparison(Metric="Alignment", Significance="***"),
            ],
        }
        unit, note, comparisons = plots._plot_statistics(
            report_data(statistics), FN_METRIC, True
        )
        self.assertEqual(unit, "image")
        self.assertIn("within a well are dependent", note)
        self.assertEqual(len(comparisons), 1)
        layout, _ = plots._comparison_layout(
            comparisons, ["Control", "Treatment"]
        )
        self.assertEqual(layout[0]["Annotation"], "ns")

    def test_nonoverlapping_blocks_reuse_bracket_levels(self):
        comparisons = [
            comparison("B", Control="A"),
            comparison("C", Control="A"),
            comparison("E", Control="D"),
            comparison("F", Control="D"),
        ]
        layout, levels = plots._comparison_layout(
            comparisons, ["A", "B", "C", "D", "E", "F"]
        )
        self.assertEqual(levels, 2)
        self.assertEqual(
            [row["Bracket_Level"] for row in layout], [0, 1, 0, 1]
        )

    def test_missing_conditions_remain_in_manifest_without_brackets(self):
        for overrides in (
            {"Control": "Unimaged control"},
            {"Treatment": "Unimaged treatment"},
        ):
            with self.subTest(overrides=overrides):
                row = comparison(
                    Status="Not tested", Significance="", **overrides
                )
                statistics = {"unit": "well", "comparisons": [row]}
                manifest, figure = self.render(report_data(statistics))
                self.assertEqual(len(figure.axes), 1)
                self.assertEqual(
                    manifest["group_order"],
                    ["Control", "Treatment", "Empty"],
                )
                self.assertEqual(figure.axes[0].get_ylim(), (0, 100))
                annotations = manifest["statistical_comparisons"]
                self.assertEqual(len(annotations), 1)
                self.assertEqual(annotations[0]["Annotation"], "Not tested")
                for key in ("Bracket_Level", "Bracket_Left", "Bracket_Right"):
                    self.assertIsNone(annotations[0][key])
                self.assertIn(
                    "no source images", annotations[0]["Annotation_Reason"]
                )
                self.assertIn("Statistics sheet", manifest["statistics_note"])

    def test_missing_condition_does_not_hide_other_comparisons(self):
        statistics = {
            "unit": "well",
            "comparisons": [
                comparison(),
                comparison(
                    "Unimaged treatment",
                    Status="Not tested",
                    Significance="",
                ),
            ],
        }
        manifest, figure = self.render(report_data(statistics))
        self.assertEqual(len(manifest["statistical_comparisons"]), 2)
        self.assertEqual(len(figure.axes), 2)
        self.assertEqual(
            [text.get_text() for text in figure.axes[0].texts], ["**"]
        )


if __name__ == "__main__":
    unittest.main()
