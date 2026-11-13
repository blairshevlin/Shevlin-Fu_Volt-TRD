"""Check the rendered figure for clipped or overlapping plot labels."""

import unittest

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.text import Text

import assemble_figure as figure_code


def descendants(artist):
    yield artist
    for child in artist.get_children():
        yield from descendants(child)


def bounds(text, renderer):
    patch = text.get_bbox_patch()
    return (patch if patch is not None else text).get_window_extent(renderer)


class FigureLayoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tables = figure_code.load_tables()
        (
            cls.figure,
            cls.axes,
            cls.right,
            cls.x_labels,
            cls.y_labels,
            cls.panel_b,
        ) = figure_code.build_figure(tables)
        cls.tables = tables
        cls.figure.canvas.draw()
        cls.renderer = cls.figure.canvas.get_renderer()

    @classmethod
    def tearDownClass(cls):
        plt.close(cls.figure)

    def assert_inside_figure(self, box):
        canvas = self.figure.bbox
        self.assertGreaterEqual(box.x0, canvas.x0)
        self.assertLessEqual(box.x1, canvas.x1)
        self.assertGreaterEqual(box.y0, canvas.y0)
        self.assertLessEqual(box.y1, canvas.y1)

    def test_panel_b_is_visible_and_above_plots(self):
        box = self.panel_b.get_window_extent(self.renderer)
        self.assert_inside_figure(box)
        self.assertGreater(box.y0, self.axes[0].bbox.y1)

    def test_x_label_chips_and_ticks_do_not_overlap(self):
        all_label_boxes = []
        for ax, label in zip(self.axes, self.x_labels):
            parts = [child for child in descendants(label) if isinstance(child, Text)]
            boxes = [bounds(part, self.renderer) for part in parts]
            all_label_boxes.append(boxes)
            ticks = [
                tick.get_window_extent(self.renderer)
                for tick in ax.get_xticklabels()
                if tick.get_visible() and tick.get_text()
            ]
            y_ticks = [
                tick.get_window_extent(self.renderer)
                for tick in ax.get_yticklabels()
                if tick.get_visible() and tick.get_text()
            ]
            for box in boxes:
                self.assert_inside_figure(box)
                for tick in ticks:
                    self.assertFalse(box.overlaps(tick), "x label overlaps a tick")
            for left, right in zip(boxes, boxes[1:]):
                self.assertFalse(left.overlaps(right), "x label segments overlap")
            for left, right in zip(ticks, ticks[1:]):
                self.assertFalse(left.overlaps(right), "x tick labels overlap")
            for tick in ticks:
                for y_tick in y_ticks:
                    self.assertFalse(tick.overlaps(y_tick), "x and y tick labels overlap")
        for left, right in zip(all_label_boxes, all_label_boxes[1:]):
            for left_box in left:
                for right_box in right:
                    self.assertFalse(left_box.overlaps(right_box), "adjacent x labels overlap")

    def test_y_labels_clear_ticks_and_plot_frames(self):
        axes_with_labels = [*self.axes, self.right]
        for ax, labels in zip(axes_with_labels, self.y_labels):
            boxes = [bounds(label, self.renderer) for label in labels]
            ticks = [
                tick.get_window_extent(self.renderer)
                for tick in ax.get_yticklabels()
                if tick.get_visible() and tick.get_text()
            ]
            for box in boxes:
                self.assert_inside_figure(box)
                for tick in ticks:
                    self.assertFalse(box.overlaps(tick), "y label overlaps a tick")
                for plot in self.axes:
                    self.assertFalse(box.overlaps(plot.bbox), "y label overlaps a plot")
            for lower, upper in zip(boxes, boxes[1:]):
                self.assertFalse(lower.overlaps(upper), "y label segments overlap")

    def test_r_squared_callouts_inside_axes_and_clear_data(self):
        axes = [*self.axes[:3], self.right]
        for index, ax in enumerate(axes):
            callouts = [text for text in ax.texts if text.get_text().startswith("R\u00b2")]
            self.assertEqual(len(callouts), 1)
            box = bounds(callouts[0], self.renderer)
            self.assert_inside_figure(box)
            self.assertGreater(box.x0, ax.bbox.x0)
            self.assertLess(box.x1, ax.bbox.x1)
            self.assertGreater(box.y0, ax.bbox.y0)
            self.assertLess(box.y1, ax.bbox.y1)
            plots = (ax,) if index < 3 else (self.axes[3], self.right)
            for plot in plots:
                for line in plot.lines:
                    if line.get_label() == "_nolegend_":
                        continue
                    if line.get_path().transformed(line.get_transform()).intersects_bbox(box):
                        self.fail("R² callout intersects plot data or fit")
                for container in plot.containers:
                    for collection in container.lines[2]:
                        for segment in collection.get_segments():
                            points = collection.get_transform().transform(segment)
                            if len(points) and box.overlaps(
                                figure_code.matplotlib.transforms.Bbox.from_extents(
                                    points[:, 0].min(), points[:, 1].min(),
                                    points[:, 0].max(), points[:, 1].max()
                                )
                            ):
                                self.fail("R² callout intersects error bars")

    def test_error_bars_match_sem_and_are_on_top(self):
        for ax, primary in zip(self.axes[:3], ("DA", "5HT", "NE")):
            table = self.tables[primary]
            x = table[f"{primary} Mean True"].to_numpy(dtype=float)
            self.assertEqual(len(ax.containers), 4)
            for container, analyte in zip(ax.containers, ("pH", "NE", "DA", "5HT")):
                prefix = analyte if analyte == primary else f"{analyte} Residual"
                y = table[f"{prefix} Mean Pred"].to_numpy(dtype=float)
                sem = table[f"{prefix} SEM"].to_numpy(dtype=float)
                self.assert_sem(container, x, y, sem, ax.transData)

        table = self.tables["pH"].dropna(subset=["pH Mean True", "pH Mean Pred"])
        x = table["pH Mean True"].to_numpy(dtype=float)
        self.assertEqual(len(self.right.containers), 4)
        for container, analyte in zip(self.right.containers, ("NE", "DA", "5HT", "pH")):
            prefix = "pH" if analyte == "pH" else f"{analyte} Residual"
            y = table[f"{prefix} Mean Pred"].to_numpy(dtype=float)
            sem = table[f"{prefix} SEM"].to_numpy(dtype=float)
            transform = self.axes[3].transData if analyte == "pH" else self.right.transData
            self.assert_sem(container, x, y, sem, transform)

    def assert_sem(self, container, x, y, sem, transform):
        marker, caps, collections = container.lines
        self.assertIsNone(marker)
        self.assertEqual(len(caps), 2)
        self.assertEqual(len(collections), 1)
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(sem)
        expected = np.stack(
            (np.column_stack((x[valid], y[valid] - sem[valid])),
             np.column_stack((x[valid], y[valid] + sem[valid]))),
            axis=1,
        )
        np.testing.assert_allclose(collections[0].get_segments(), expected)
        self.assertEqual(collections[0].get_zorder(), 6)
        self.assertEqual(collections[0].get_transform(), transform)
        self.assertTrue(all(cap.get_zorder() == 6 for cap in caps))


if __name__ == "__main__":
    unittest.main()
