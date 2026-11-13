"""Rebuild panel b of Extended Data Fig. 4 as an editable SVG and a 300-dpi TIFF.

Adapted for this repository from the authors' two-panel script: panel a (BioRender
artwork) is not rebuilt here, and the binned data are read from the figure's
source-data workbook instead of four CSVs.
"""

from pathlib import Path
import logging
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import numpy as np
import pandas as pd
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.ticker import FuncFormatter
from PIL import Image

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


ROOT = Path(__file__).resolve().parents[3]
SOURCE_DATA = ROOT / "data" / "figures" / "figureExtended4_source_data.xlsx"
SOURCE_SHEET = "Extended Data Figure 4b"
OUTPUT_SVG = ROOT / "results" / "from_source_data" / "Extended-Data_Figure4b.svg"
OUTPUT_TIFF = ROOT / "results" / "from_source_data" / "Extended-Data_Figure4b.tif"
OUTPUT_DPI = 300
PLOT_FONT_SIZE_PT = 7.5
PANEL_LETTER_SIZE_PT = 10.7
MARKER_SIZE_PT = 4.6
# Panel b's area of the published 8.5 x 4.65 in figure (panel a took the top 2.35 in)
FIGURE_SIZE_INCHES = (8.5, 2.3)

COLORS = {
    "DA": (0.796, 0.094, 0.114),
    "5HT": (0.129, 0.443, 0.710),
    "NE": (0.137, 0.545, 0.271),
    "pH": (0.533, 0.337, 0.655),
}
ANALYTES = ("DA", "5HT", "NE", "pH")

try:
    findfont(FontProperties(family="Arial"), fallback_to_default=False)
    FONT_FAMILY = "Arial"
except ValueError:
    # Matplotlib's persistent font cache can predate a new font installation.
    for font in font_manager.FontManager().ttflist:
        if font.name.casefold() == "arial":
            font_manager.fontManager.addfont(font.fname)
    try:
        findfont(FontProperties(family="Arial"), fallback_to_default=False)
        FONT_FAMILY = "Arial"
    except ValueError:
        FONT_FAMILY = "Liberation Sans"

matplotlib.rcParams.update(
    {
        "font.family": FONT_FAMILY,
        "font.sans-serif": [FONT_FAMILY, "DejaVu Sans"],
        "font.size": PLOT_FONT_SIZE_PT,
        "svg.fonttype": "none",
        "svg.image_inline": True,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 1.8,
        "ytick.major.size": 1.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "savefig.facecolor": "white",
    }
)


def load_tables():
    """Rebuilds the four per-analyte bin tables (original column names) from the long-format sheet."""
    if not SOURCE_DATA.is_file():
        raise FileNotFoundError(f"Required source data not found: {SOURCE_DATA}")
    sheet = pd.read_excel(SOURCE_DATA, sheet_name=SOURCE_SHEET)
    tables = {}
    for analyte in ANALYTES:
        rows = sheet[sheet["plot"] == analyte]
        if rows.empty:
            raise ValueError(f"No rows for {analyte} in sheet {SOURCE_SHEET!r}")
        table = rows.groupby("bin")[["true_mean", "true_sem"]].first()
        table.columns = [f"{analyte} Mean True", f"{analyte} x_SEM"]
        for other, values in rows.groupby("analyte"):
            prefix = other if other == analyte else f"{other} Residual"
            values = values.set_index("bin")
            table[f"{prefix} Mean Pred"] = values["mean"]
            table[f"{prefix} SEM"] = values["sem"]
        tables[analyte] = table
    return tables


def tint(color):
    return tuple(0.6 + 0.4 * channel for channel in color)


def text_area(text, *, color="black", font_size=PLOT_FONT_SIZE_PT, bbox=None):
    props = {
        "color": color,
        "fontsize": font_size,
        "fontfamily": FONT_FAMILY,
        "math_fontfamily": "dejavusans",
    }
    if bbox is not None:
        props["bbox"] = bbox
    return TextArea(text, textprops=props)


def analyte_chip(analyte):
    color = COLORS[analyte]
    return text_area(
        analyte,
        color="black",
        bbox={
            "boxstyle": "round,pad=0.12",
            "facecolor": tint(color),
            "edgecolor": "black",
            "linewidth": 0.45,
        },
    )


def add_x_label(ax):
    parts = [text_area("True [")]
    for index, analyte in enumerate(ANALYTES):
        if index:
            parts.append(text_area(" or "))
        parts.append(analyte_chip(analyte))
    parts.append(text_area("] (nM)"))

    label = AnchoredOffsetbox(
        loc="lower center",
        child=HPacker(children=parts, align="center", pad=0, sep=1.5),
        pad=0,
        borderpad=0,
        frameon=False,
        bbox_to_anchor=(0.5, -0.22),
        bbox_transform=ax.transAxes,
    )
    label.set_clip_on(False)
    ax.add_artist(label)
    return label


def add_y_label(ax, analytes, x_position, prefix="Predicted ["):
    font = FontProperties(family=FONT_FAMILY, size=PLOT_FONT_SIZE_PT)
    segments = [("text", prefix, None)]
    for index, analyte in enumerate(analytes):
        if index:
            segments.append(("text", " or ", None))
        segments.append(("chip", analyte, analyte))
    segments.append(("text", "] (nM)", None))

    renderer = ax.figure.canvas.get_renderer()
    measured = []
    for kind, text, analyte in segments:
        width = renderer.get_text_width_height_descent(text, font, False)[0]
        width *= 72 / ax.figure.dpi
        if kind == "chip":
            width += PLOT_FONT_SIZE_PT * 0.24
        measured.append((kind, text, analyte, width))

    axes_height_pt = ax.get_position().height * ax.figure.get_figheight() * 72
    gap_pt = PLOT_FONT_SIZE_PT * 0.24
    total_width_pt = sum(item[-1] for item in measured) + gap_pt * (len(measured) - 1)
    cursor = 0.5 - total_width_pt / (2 * axes_height_pt)

    labels = []
    for kind, text, analyte, width in measured:
        center = cursor + width / (2 * axes_height_pt)
        props = {"rotation": 90, "ha": "center", "va": "center", "clip_on": False}
        if kind == "chip":
            color = COLORS[analyte]
            label = ax.text(
                x_position,
                center,
                text,
                transform=ax.transAxes,
                fontsize=PLOT_FONT_SIZE_PT,
                color="black",
                bbox={
                    "boxstyle": "round,pad=0.12",
                    "facecolor": tint(color),
                    "edgecolor": "black",
                    "linewidth": 0.45,
                },
                **props,
            )
        else:
            label = ax.text(
                x_position,
                center,
                text,
                transform=ax.transAxes,
                fontsize=PLOT_FONT_SIZE_PT,
                color="black",
                **props,
            )
        labels.append(label)
        cursor += (width + gap_pt) / axes_height_pt
    return labels


def add_markers(ax, x, y, analyte, zorder):
    ax.plot(
        x,
        y,
        linestyle="none",
        marker="o",
        markerfacecolor=COLORS[analyte],
        markeredgecolor="black",
        markeredgewidth=0.55,
        markersize=MARKER_SIZE_PT,
        zorder=zorder,
    )


def add_sem_bars(ax, x, y, sem, *, transform=None):
    ax.errorbar(
        x,
        y,
        yerr=sem,
        fmt="none",
        ecolor="black",
        elinewidth=0.6,
        capsize=2,
        capthick=0.6,
        zorder=6,
        transform=transform or ax.transData,
    )


def plot_analyte(ax, table, analyte):
    true_column = f"{analyte} Mean True"
    x = table[true_column].to_numpy(dtype=float)

    series = []
    for other in ("pH", "NE", "DA", "5HT"):
        is_primary = other == analyte
        value_column = (
            f"{other} Mean Pred"
            if is_primary
            else f"{other} Residual Mean Pred"
        )
        sem_column = f"{other} SEM" if is_primary else f"{other} Residual SEM"
        if value_column not in table:
            continue
        y = table[value_column].to_numpy(dtype=float)
        sem = table[sem_column].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(sem)
        if not valid.any():
            continue
        series.append((other, x[valid], y[valid], sem[valid], is_primary))

    primary = next(item for item in series if item[-1])
    _, fit_x, fit_y, _, _ = primary
    slope, intercept = np.polyfit(fit_x, fit_y, 1)
    fitted = slope * fit_x + intercept
    residual_sum_squares = np.sum((fit_y - fitted) ** 2)
    total_sum_squares = np.sum((fit_y - np.mean(fit_y)) ** 2)
    r_squared = 1 - residual_sum_squares / total_sum_squares

    # Keep the primary series on top so its colored markers remain visible.
    for other, values_x, values_y, sem, is_primary in series:
        add_markers(ax, values_x, values_y, other, 3 if is_primary else 2)
    ax.plot(fit_x, fitted, color="black", linewidth=0.75, zorder=4)
    for _, values_x, values_y, sem, _ in series:
        add_sem_bars(ax, values_x, values_y, sem)

    ax.set_xlim((-300, 2800))
    ax.set_ylim((-800, 2800))
    ax.set_xticks((0, 1000, 2000))
    ax.set_xticks((500, 1500, 2500), minor=True)
    ax.tick_params(axis="x", which="minor", direction="in", length=1.8, width=0.6)
    ax.set_yticks(np.arange(0, 2501, 500))
    ax.tick_params(labelsize=PLOT_FONT_SIZE_PT, pad=2)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color("black")

    label = f"R\u00b2 = {np.floor(r_squared * 1000) / 1000:.3f}"
    ax.text(
        0.08,
        0.94,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=PLOT_FONT_SIZE_PT,
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": tint(COLORS[analyte]),
            "edgecolor": "black",
            "linewidth": 0.65,
        },
        clip_on=False,
        zorder=7,
    )


def plot_ph(ax, table):
    table = table.dropna(subset=["pH Mean True", "pH Mean Pred"])
    x = table["pH Mean True"].to_numpy(dtype=float)
    y = table["pH Mean Pred"].to_numpy(dtype=float)
    sem = table["pH SEM"].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    residual_sum_squares = np.sum((y - fitted) ** 2)
    total_sum_squares = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - residual_sum_squares / total_sum_squares

    add_markers(ax, x, y, "pH", 3)
    ax.plot(x, fitted, color="black", linewidth=0.75, zorder=4)

    ax.set_xlim((6.5, 8.2))
    ax.set_ylim((6.5, 8.2))
    ax.set_xticks((6.5, 7, 7.5, 8))
    ax.set_yticks((6.5, 7, 7.5, 8))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax.tick_params(labelsize=PLOT_FONT_SIZE_PT, pad=2)
    ax.tick_params(axis="x", pad=5)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color("black")

    right = ax.twinx()
    right.set_ylim((-800, 2800))
    right.set_yticks(np.arange(0, 2501, 500))
    right.tick_params(labelsize=PLOT_FONT_SIZE_PT, pad=2)
    residual_series = []
    for other in ("NE", "DA", "5HT"):
        residual = table[f"{other} Residual Mean Pred"].to_numpy(dtype=float)
        residual_sem = table[f"{other} Residual SEM"].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(residual) & np.isfinite(residual_sem)
        residual_series.append((x[valid], residual[valid], residual_sem[valid]))
        add_markers(right, x[valid], residual[valid], other, 2)
    for values_x, values_y, values_sem in residual_series:
        add_sem_bars(right, values_x, values_y, values_sem)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(sem)
    add_sem_bars(right, x[valid], y[valid], sem[valid], transform=ax.transData)
    for spine in right.spines.values():
        spine.set_color("black")
    right.set_facecolor("none")
    right.patch.set_visible(False)

    label = f"R\u00b2 = {np.floor(r_squared * 1000) / 1000:.3f}"
    right.text(
        0.08,
        0.94,
        label,
        transform=right.transAxes,
        ha="left",
        va="top",
        fontsize=PLOT_FONT_SIZE_PT,
        bbox={
            "boxstyle": "round,pad=0.22",
            "facecolor": tint(COLORS["pH"]),
            "edgecolor": "black",
            "linewidth": 0.65,
        },
        clip_on=False,
        zorder=7,
    )
    return right


def build_figure(tables):
    figure = plt.figure(figsize=FIGURE_SIZE_INCHES, dpi=OUTPUT_DPI, facecolor="white")

    plot_bottom = 0.43 / FIGURE_SIZE_INCHES[1]
    plot_width = 0.155
    plot_height = plot_width * FIGURE_SIZE_INCHES[0] / FIGURE_SIZE_INCHES[1]
    plot_lefts = (0.065, 0.29, 0.515, 0.74)
    axes = [
        figure.add_axes((left, plot_bottom, plot_width, plot_height))
        for left in plot_lefts
    ]

    for ax, analyte in zip(axes[:3], ("DA", "5HT", "NE")):
        plot_analyte(ax, tables[analyte], analyte)
    ph_right = plot_ph(axes[3], tables["pH"])

    y_labels = [
        add_y_label(ax, (analyte,), x_position=-0.29)
        for ax, analyte in zip(axes, ANALYTES)
    ]
    y_labels.append(add_y_label(ph_right, ("DA", "5HT", "NE"), x_position=1.32))
    x_labels = [add_x_label(ax) for ax in axes]
    panel_b_label = figure.text(
        0.012,
        2.16 / FIGURE_SIZE_INCHES[1],
        "b",
        fontsize=PANEL_LETTER_SIZE_PT,
        ha="left",
        va="top",
    )
    return figure, axes, ph_right, x_labels, y_labels, panel_b_label


def main():
    tables = load_tables()
    figure, _, _, _, _, _ = build_figure(tables)

    figure.savefig(OUTPUT_TIFF, dpi=OUTPUT_DPI, pil_kwargs={"compression": "tiff_lzw"})
    figure.savefig(OUTPUT_SVG, format="svg", dpi=OUTPUT_DPI)
    plt.close(figure)
    svg = OUTPUT_SVG.read_text(encoding="utf-8")
    OUTPUT_SVG.write_text(
        svg.replace("Liberation Sans", "Arial"),
        encoding="utf-8",
    )
    with Image.open(OUTPUT_TIFF) as image:
        rgb = image.convert("RGB")
        rgb.save(
            OUTPUT_TIFF,
            format="TIFF",
            dpi=(OUTPUT_DPI, OUTPUT_DPI),
            compression="tiff_lzw",
        )
    print(f"Saved {OUTPUT_TIFF} and editable {OUTPUT_SVG} at {OUTPUT_DPI} dpi")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Figure assembly failed: {error}", file=sys.stderr)
        raise
