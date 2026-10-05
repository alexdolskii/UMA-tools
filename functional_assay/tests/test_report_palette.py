"""Render roles, missing days, hatches and unchanged measurement values."""

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from functional_assay import plot_palette, report_output, survival_output
from functional_assay.survival_data import calculate_changes, observations
from matplotlib.colors import to_rgba
from matplotlib.markers import MarkerStyle


def synthetic_data(count=6, control=True):
    """A fabricated plate with three technical wells per condition."""
    groups = ["Control"] + [f"Treatment {i}" for i in range(1, count)]
    design, rows = [], []
    for group_index, group in enumerate(groups):
        for index in range(3):
            position = group_index * 3 + index
            annotation = {
                "Well": f"{'ABCDEFGH'[position // 12]}{position % 12 + 1:02}",
                "Comparison_Block": "Block_01",
                "Group": group,
                "Is_Control": control and group_index == 0,
            }
            design.append(annotation)
            for day in (1, 3, 5):
                value = 10 + group_index * 7 + index * 3 + day * (index + 1)
                rows.append(
                    {
                        **annotation,
                        "Day": day,
                        "Object_Count": value,
                        "Mask_Area_um2": value * 12.5,
                        "Annotation_Status": "MAPPED",
                    }
                )
    data = {
        "blocks": [
            {
                "id": "Block_01",
                "groups": groups,
                "control": "Control" if control else None,
                "rgb": "#ABCDEF",
            }
        ],
        "design": design,
        "rows": rows,
        "days": [1, 3, 5],
        "baseline_day": 1,
        "difference_days": [3, 5],
        "comparisons": [],
        "statistics_enabled": False,
    }
    data["changes"] = calculate_changes(rows, data, 1, [3, 5])
    return data


class ReportPaletteTests(unittest.TestCase):
    def test_six_condition_boundary_uses_full_template(self):
        for count in (2, 5, 6, 8, 9):
            with self.subTest(conditions=count):
                data = synthetic_data(count)
                before = plot_palette.prepare_palette(data, survival=True)
                styles = before["conditions"]["Block_01"]
                self.assertEqual(styles["Control"]["color"], "#B5B1D8")
                self.assertEqual(styles["Treatment 1"]["color"], "#004F46")
                self.assertEqual(
                    any(s["hatch"] for s in styles.values()), count >= 6
                )
                if count >= 6:
                    self.assertEqual(
                        len({s["hatch"] for s in styles.values()}), count
                    )
                self.assertEqual(
                    styles["Control"]["palette_mode"],
                    "historical" if count <= 8 else "green_tints",
                )
                # An absent condition and a wholly missing day cannot
                # reassign any color, hatch, well marker or jitter slot.
                data["rows"] = [
                    r
                    for r in data["rows"]
                    if r["Group"] != "Control" and r["Day"] != 3
                ]
                self.assertEqual(
                    plot_palette.prepare_palette(data, survival=True), before
                )
                ordinary = plot_palette.prepare_palette(data)
                self.assertTrue(
                    all(
                        not s["hatch"]
                        for s in ordinary["conditions"]["Block_01"].values()
                    )
                )

    def test_lavender_without_bold_never_creates_a_control(self):
        for count in (2, 8, 9):
            data = synthetic_data(count, control=False)
            styles = plot_palette.prepare_palette(data)["conditions"][
                "Block_01"
            ]
            self.assertFalse(any(s["is_control"] for s in styles.values()))
            self.assertEqual(styles["Control"]["color"], "#004F46")
            if count <= 8:
                self.assertEqual(styles["Treatment 1"]["color"], "#B5B1D8")
            else:
                self.assertTrue(
                    all(
                        s["palette_mode"] == "green_tints"
                        and s["color"] != "#B5B1D8"
                        for s in styles.values()
                    )
                )
            self.assertIsNone(data["blocks"][0]["control"])

    def test_repeated_condition_name_can_have_different_roles_in_two_blocks(
        self,
    ):
        data = synthetic_data(2)
        data["blocks"].append(
            {
                "id": "Block_02",
                "groups": ["Control", "Reference"],
                "control": "Reference",
                "rgb": "#FEDCBA",
            }
        )
        for well, group in [("H01", "Control"), ("H02", "Reference")]:
            data["design"].append(
                {
                    "Well": well,
                    "Group": group,
                    "Comparison_Block": "Block_02",
                    "Is_Control": group == "Reference",
                }
            )
        palette = plot_palette.prepare_palette(data)
        styles = palette["conditions"]
        self.assertEqual(styles["Block_01"]["Control"]["color"], "#B5B1D8")
        self.assertEqual(styles["Block_02"]["Control"]["color"], "#004F46")
        self.assertEqual(styles["Block_02"]["Reference"]["color"], "#B5B1D8")
        self.assertEqual(len(palette["wells"]), 8)

    def assert_point(self, collection, value, marker):
        np.testing.assert_array_equal(
            collection.get_facecolors()[0], to_rgba("#D0D0D0")
        )
        np.testing.assert_array_equal(
            collection.get_edgecolors()[0], to_rgba("#111314")
        )
        self.assertEqual(collection.get_offsets()[0, 1], value)
        shape = MarkerStyle(marker)
        np.testing.assert_array_equal(
            collection.get_paths()[0].vertices,
            shape.get_path().transformed(shape.get_transform()).vertices,
        )

    def test_six_survival_figures_keep_values_and_missing_well_shapes(self):
        data = synthetic_data(6)
        data["rows"] = [
            row
            for row in data["rows"]
            if not (row["Day"] == 3 and row["Well"] == "A01")
            and row["Day"] != 5
        ]
        data["changes"] = calculate_changes(data["rows"], data, 1, [3, 5])
        before = copy.deepcopy(data)
        figures = []
        with tempfile.TemporaryDirectory() as folder:
            with patch(
                "matplotlib.figure.Figure.savefig",
                lambda figure, *a, **k: figures.append(figure),
            ):
                plots = survival_output.render_plots(
                    data, Path(folder), "Test"
                )
        self.assertEqual(len(plots), 6)
        self.assertEqual(data["rows"], before["rows"])
        self.assertEqual(data["changes"], before["changes"])
        for plot, figure in zip(plots, figures):
            axes = [a for a in figure.axes if a.collections]
            self.assertEqual(len(axes), 1)
            axis = axes[0]
            source = {
                (r["Day"], r["Well"]): r["Value"]
                for r in observations(before, plot["view"], plot["metric"])
            }
            self.assertEqual(len(axis.collections), len(source))
            for point, collection in zip(plot["points"], axis.collections):
                key = (point["Day"], point["Well"])
                index = (
                    next(
                        i
                        for i, r in enumerate(data["design"])
                        if r["Well"] == point["Well"]
                    )
                    % 3
                )
                self.assert_point(
                    collection, source[key], ("o", "s", "^")[index]
                )
                # The second well retains its central slot even when
                # the first well or a complete day is absent.
                if index == 1:
                    x = collection.get_offsets()[0, 0]
                    self.assertEqual(x, round(x))
            # Exact colors and pattern sequence, not transparent Excel fills.
            colors = [
                "#B5B1D8",
                "#004F46",
                "#F37420",
                "#051230",
                "#80719E",
                "#EBD3A2",
            ]
            hatches = ["", "//", "xx", "..", "\\\\", "++"]
            for index, box in enumerate(axis.patches):
                np.testing.assert_array_equal(
                    box.get_facecolor(), to_rgba(colors[index % 6])
                )
                self.assertEqual(box.get_hatch(), hatches[index % 6])

    def test_two_single_day_plots_use_ordinary_lavender_without_hatching(self):
        data = synthetic_data(2, control=False)
        data["rows"] = [r for r in data["rows"] if r["Day"] == 1]
        figures = []
        with tempfile.TemporaryDirectory() as folder:
            with patch(
                "matplotlib.figure.Figure.savefig",
                lambda figure, *a, **k: figures.append(figure),
            ):
                plots = report_output.render_plots(data, Path(folder), "Test")
        self.assertEqual(len(plots), 2)
        for plot, figure in zip(plots, figures):
            axis = next(a for a in figure.axes if a.collections)
            for color, box in zip(["#004F46", "#B5B1D8"], axis.patches):
                np.testing.assert_array_equal(
                    box.get_facecolor(), to_rgba(color)
                )
                self.assertFalse(box.get_hatch())
            for row, collection in zip(data["rows"], axis.collections):
                index = int(row["Well"][1:]) - 1
                self.assert_point(
                    collection, row[plot["metric"]], ("o", "s", "^")[index % 3]
                )


if __name__ == "__main__":
    unittest.main()
