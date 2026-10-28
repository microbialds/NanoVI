"""Shared figure style for the NanoVI manuscript (decision D12, journal-neutral).

Sizes, fonts, line widths, the Okabe-Ito palette, the fixed colour of every
tool/database, the panel-label helper and the writer that saves each figure as
PDF, SVG, PNG (600 dpi) and its plotted data as TSV.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

MM = 1 / 25.4
SINGLE = 85 * MM   # single column, inches
DOUBLE = 175 * MM  # double column, inches
MAX_HEIGHT = 200 * MM

FONT_SIZE = 7
PANEL_LABEL_SIZE = 8
LINE_WIDTH = 0.6
DPI = 600

OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}

# One colour per tool and database, the same in every figure.
# Keys are the (tool, database) values of fig_comparators.tsv.
TOOLS = {
    ("nanovi", "nanovi_exact"): ("NanoVI (GTDB r226)", OKABE_ITO["blue"]),
    ("nanovi", "emu_default"): ("NanoVI (Emu db)", OKABE_ITO["sky_blue"]),
    ("emu", "emu_default"): ("Emu (Emu db)", OKABE_ITO["vermillion"]),
    ("emu", "nanovi_exact"): ("Emu (GTDB r226)", OKABE_ITO["orange"]),
    ("wf16s_minimap2", "ncbi_16s_18s"): ("wf-16s minimap2 (NCBI)", OKABE_ITO["bluish_green"]),
    ("wf16s_kraken2", "ncbi_16s_18s"): ("wf-16s kraken2 (NCBI)", OKABE_ITO["reddish_purple"]),
    ("kraken2_bracken", "ncbi_16s_18s"): ("Kraken2 + Bracken (NCBI)", OKABE_ITO["black"]),
    ("nanoclust", "ncbi_16S_ribosomal_RNA"): ("NanoCLUST (NCBI type strains)", OKABE_ITO["yellow"]),
}
NANOVI = TOOLS[("nanovi", "nanovi_exact")][1]
EMU_OWN_DB = TOOLS[("emu", "emu_default")][1]
EMU_GTDB = TOOLS[("emu", "nanovi_exact")][1]
# Light colours get a dark edge so that they stay visible on white.
EDGE_FOR_LIGHT = "#333333"
LIGHT_COLOURS = {OKABE_ITO["yellow"]}


def edge_for(colour):
    return EDGE_FOR_LIGHT if colour in LIGHT_COLOURS else colour


def mix(colour, other, t):
    """Blend colour towards other by t (0 = colour, 1 = other)."""
    a, b = to_rgb(colour), to_rgb(other)
    return tuple((1 - t) * x + t * y for x, y in zip(a, b))


# Read depth (simulations): one hue, the NanoVI blue, from light to dark.
DEPTH_COLOURS = {
    1000: mix(NANOVI, "white", 0.55),
    10000: NANOVI,
    50000: mix(NANOVI, "black", 0.55),
}
DEPTH_MARKERS = {1000: "o", 10000: "s", 50000: "D"}

# Simulated samples: marker by community shape, filled for R10.4.1, open for R9.4.1.
SHAPE_MARKERS = {"even": "o", "log": "^"}
CHEMISTRY_LABELS = {"r1041": "R10.4.1", "r941": "R9.4.1"}
SHAPE_LABELS = {"even": "even", "log": "log-staggered"}


def pick_font():
    """Arial when installed, DejaVu Sans otherwise (D12)."""
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
        return "Arial"
    except ValueError:
        return "DejaVu Sans"


def apply_style():
    font = pick_font()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "DejaVu Sans"],
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
        "legend.title_fontsize": FONT_SIZE,
        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.3,
        "axes.linewidth": LINE_WIDTH,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelpad": 2.5,
        "xtick.major.width": LINE_WIDTH,
        "ytick.major.width": LINE_WIDTH,
        "xtick.minor.width": LINE_WIDTH * 0.8,
        "ytick.minor.width": LINE_WIDTH * 0.8,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.major.pad": 2,
        "ytick.major.pad": 2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "patch.linewidth": 0.5,
        "hatch.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "nanovi",
        "savefig.dpi": DPI,
        "figure.dpi": 150,
        "mathtext.default": "regular",
    })
    return font


def figure(width, height_mm):
    height = min(height_mm * MM, MAX_HEIGHT)
    return plt.figure(figsize=(width, height), layout="constrained")


def panel_labels(fig, pairs):
    """Bold lower-case panel labels at the top-left of each panel.

    Placed after a first layout pass, at the left edge of the panel's tick labels
    and y label, just above the axes.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, letter in pairs:
        tight = ax.get_tightbbox(renderer)
        box = ax.get_window_extent(renderer)
        x = (tight.x0 - box.x0) / box.width
        y = max(1.0, (tight.y1 - box.y0) / box.height)
        ax.text(x, y, letter, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=PANEL_LABEL_SIZE, fontweight="bold")


def italic_ticklabels(labels):
    for t in labels:
        t.set_fontstyle("italic")


def save(fig, name, data, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext, meta in (("pdf", {"CreationDate": None, "ModDate": None}),
                      ("svg", {"Date": None}),
                      ("png", {})):
        path = out_dir / f"{name}.{ext}"
        fig.savefig(path, metadata=meta or None, dpi=DPI if ext == "png" else None)
        written.append(path)
    path = out_dir / f"{name}_data.tsv"
    data.to_csv(path, sep="\t", index=False)
    written.append(path)
    plt.close(fig)
    return written
