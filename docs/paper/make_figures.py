#!/usr/bin/env python3
"""Draw the NanoVI manuscript figures from the benchmark tables.

Every plotted value is read from a table in docs/benchmarks/data/ (columns are
described in docs/benchmarks/methods_2026.md). The only constants drawn are
reference marks: the nominal 0.95 coverage, zero error and the y = x diagonal.

    python docs/paper/make_figures.py --out internal_docs/manuscript/figures
    python docs/paper/make_figures.py --out internal_docs/manuscript/figures --only fig3

For each figure it writes <name>.pdf, <name>.svg, <name>.png (600 dpi) and
<name>_data.tsv (the plotted rows, with the source file of each row).
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402
from style import plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "docs" / "benchmarks" / "data"


def read(name):
    return pd.read_csv(DATA / name, sep="\t")


def with_source(df, name, **cols):
    df = df.copy()
    for k, v in cols.items():
        df[k] = v
    df["source_file"] = name
    return df


def depth_label(depth):
    depth = int(depth)
    return f"{depth // 1000}k" if depth >= 1000 and depth % 1000 == 0 else str(depth)


def sample_marker(chemistry, shape, colour):
    """Marker kwargs: shape gives the symbol, R10.4.1 filled, R9.4.1 open."""
    filled = chemistry == "r1041"
    return dict(marker=style.SHAPE_MARKERS[shape], markerfacecolor=colour if filled else "white",
                markeredgecolor=colour, markeredgewidth=0.7, linestyle="none")


def short_species(name):
    genus, _, rest = name.partition(" ")
    return f"{genus[0]}. {rest}" if rest else name


def fmt_number(v):
    v = float(v)
    return f"{int(v):,}" if v.is_integer() else f"{v:,g}"


# ─── Figure 1: database ─────────────────────────────────────────────────────

def fig1_database():
    hist_name = "db_audit_ssu_all_seqs_per_species_hist.tsv"
    summ_name = "db_audit_ssu_all_summary.tsv"
    dbn_name = "fig_db_N.tsv"
    hist = read(hist_name)
    summ = read(summ_name).set_index("metric")["value"]
    marks = [("median", "sequences per species: median"),
             ("95th percentile", "sequences per species: 95th percentile"),
             ("maximum", "sequences per species: maximum")]
    marks = [(label, metric, float(summ[metric])) for label, metric in marks]

    # log2 bins: [1], [2, 3], [4, 7], ...
    b = np.floor(np.log2(hist["db_sequences"])).astype(int)
    binned = hist.groupby(b)["species"].sum()
    lo = 2.0 ** binned.index.values
    hi = 2.0 ** (binned.index.values + 1) - 1

    fig = style.figure(style.DOUBLE, 62)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax.bar(lo, binned.values, width=(hi + 1 - lo) * 0.86, align="edge", color=style.NANOVI,
           edgecolor="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.7, (hi.max() + 1) * 1.6)
    ax.set_ylim(0.6, binned.max() * 15)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_xlabel("Database sequences per species")
    ax.set_ylabel("Species")
    ys = [0.97, 0.87, 0.97]
    for (label, _, value), y in zip(marks, ys):
        ax.axvline(value, color="#555555", lw=0.6, ls="--", zorder=0)
        right = value < hi.max() / 4
        ax.text(value * (1.12 if right else 1 / 1.12), y, f"{label}\n{fmt_number(value)}",
                transform=ax.get_xaxis_transform(), ha="left" if right else "right", va="top",
                fontsize=style.FONT_SIZE)
    data_a = pd.DataFrame({"panel": "a", "bin_low": lo.astype(int), "bin_high": hi.astype(int),
                           "species": binned.values, "source_file": hist_name})
    data_a = pd.concat([data_a, pd.DataFrame({"panel": "a", "metric": [m for _, m, _ in marks],
                                              "value": [v for _, _, v in marks],
                                              "source_file": summ_name})])

    # (b) species and genus L1 per database variant and -N
    dbn = read(dbn_name)
    arms = [("exact", 3), ("exact", 10), ("exact", 50), ("c995", 3), ("reps", 3)]
    sel = pd.concat([dbn[(dbn["db"] == db) & (dbn["N"] == n)] for db, n in arms])
    samples = list(dict.fromkeys(sel["sample"]))
    offsets = np.linspace(-0.18, 0.18, len(samples))
    axb = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
    ymax = sel[["species_l1", "genus_l1"]].max().max() * 1.08
    for a, metric, ylabel in zip(axb, ["species_l1", "genus_l1"],
                                 ["Species L1 distance to truth", "Genus L1 distance to truth"]):
        for i, (db, n) in enumerate(arms):
            d = sel[(sel["db"] == db) & (sel["N"] == n)]
            for j, s in enumerate(samples):
                row = d[d["sample"] == s]
                shape, chem = s.split("_")[:2]
                a.plot(i + offsets[j], row[metric].values, ms=3.6,
                       **sample_marker(chem, shape, style.NANOVI))
            a.plot([i - 0.3, i + 0.3], [d[metric].mean()] * 2, color="black", lw=1.0)
        a.set_xticks(range(len(arms)), [f"{db}\n-N {n}" for db, n in arms])
        a.set_xlim(-0.6, len(arms) - 0.4)
        a.set_ylim(0, ymax)
        a.set_ylabel(ylabel)
    handles = []
    for s in samples:
        shape, chem = s.split("_")[:2]
        handles.append(Line2D([], [], ms=3.6, **sample_marker(chem, shape, style.NANOVI),
                              label=f"{style.SHAPE_LABELS[shape]}, {style.CHEMISTRY_LABELS[chem]}"))
    handles.append(Line2D([], [], color="black", lw=1.0, label="mean of the four samples"))
    axb[1].legend(handles=handles, loc="upper left", handletextpad=0.3)
    style.panel_labels(fig, [(ax, "a"), (axb[0], "b")])

    data_b = sel[["sample", "db", "N", "species_l1", "genus_l1"]].copy()
    data_b.insert(0, "panel", "b")
    data_b["source_file"] = dbn_name
    return fig, pd.concat([data_a, data_b], ignore_index=True)


# ─── Figure 2: comparators ──────────────────────────────────────────────────

def fig2_comparators():
    name = "fig_comparators.tsv"
    c = read(name)
    real = c[c["kind"] == "real"]
    keys = [k for k in style.TOOLS if ((real["tool"] == k[0]) & (real["database"] == k[1])).any()]
    groups = {k: real[(real["tool"] == k[0]) & (real["database"] == k[1])] for k in keys}
    order = sorted(keys, key=lambda k: groups[k]["l1_group"].mean())
    n_common = pd.Series([len(groups[k]) for k in keys]).mode().iloc[0]
    rng = np.random.default_rng(2)

    fig = style.figure(style.DOUBLE, 78)
    gs = fig.add_gridspec(1, 2)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1], sharey=axa)
    rows = []
    for i, k in enumerate(order):
        y = len(order) - 1 - i
        label, colour = style.TOOLS[k]
        g = groups[k].copy()
        g["mock"] = np.where(g["sample"].str.lower().str.contains("gutmock"), "gut", "zymo")
        jit = rng.uniform(-0.2, 0.2, len(g))
        for ax, metric in ((axa, "l1_group"), (axb, "group_fp_count")):
            for mock, marker in (("zymo", "o"), ("gut", "s")):
                m = (g["mock"] == mock).values
                ax.scatter(g.loc[m, metric], y + jit[m], s=11, marker=marker, color=colour,
                           edgecolors=style.edge_for(colour), linewidths=0.4, alpha=0.85, zorder=2)
            mean = g[metric].mean()
            ax.plot([mean, mean], [y - 0.38, y + 0.38], color="black", lw=1.1, zorder=3,
                    path_effects=[pe.Stroke(linewidth=2.6, foreground="white"), pe.Normal()])
        g["label"] = label
        g["y_order"] = i + 1
        rows.append(g)
    labels = []
    for k in order:
        label = style.TOOLS[k][0]
        n = len(groups[k])
        labels.append(label if n == n_common else f"{label}\n{n} samples")
    axa.set_yticks(range(len(order) - 1, -1, -1), labels)
    axa.set_ylim(-0.6, len(order) - 0.4)
    axa.set_xlim(left=0)
    axa.set_xlabel("Identical-16S group L1 distance to truth")
    axb.set_xscale("symlog", linthresh=1, linscale=0.6)
    axb.set_xlim(-0.3, real["group_fp_count"].max() * 1.6)
    ticks = [t for t in (0, 1, 2, 5, 10, 20, 50, 100) if t <= real["group_fp_count"].max() * 1.6]
    axb.set_xticks(ticks, [str(t) for t in ticks])
    axb.xaxis.set_minor_locator(plt.NullLocator())
    axb.set_xlabel("False-positive groups per sample")
    plt.setp(axb.get_yticklabels(), visible=False)
    for ax in (axa, axb):
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color="#E5E5E5", lw=0.4, zorder=0)
        ax.set_axisbelow(True)
    handles = [Line2D([], [], marker="o", ls="none", color="#666666", ms=3.6,
                      label="Zymo community standard"),
               Line2D([], [], marker="s", ls="none", color="#666666", ms=3.6, label="gut mock"),
               Line2D([], [], color="black", lw=1.1, label="mean over samples")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3)
    style.panel_labels(fig, [(axa, "a"), (axb, "b")])

    data = pd.concat(rows)[["y_order", "label", "tool", "database", "group", "sample", "mock",
                            "l1_group", "group_fp_count"]]
    data = data.rename(columns={"group": "sample_group"})
    data["source_file"] = name
    return fig, data


# ─── Figure 3: Zymo composition ─────────────────────────────────────────────

ZYMO_SAMPLES = [
    ("preprint_zymo_SRR23636353", "SRR23636353 (preprint mock)"),
    ("emu_zymo_SRR10391201", "SRR10391201 (Emu paper mock)"),
]


def fig3_zymo_composition():
    name = "real_species.tsv"
    s = read(name)
    ref = s[s["arm"] == "reference"]
    msplus = sorted(x for x in ref["sample"].unique() if re.match(r"zymo_MSPlus_sup_rep\d+$", x))
    samples = ZYMO_SAMPLES + [(x, f"MSPlus sup rep {x.rsplit('rep', 1)[1]}") for x in msplus]
    ref = ref[ref["sample"].isin([x for x, _ in samples])]

    exp = ref[ref["in_truth"] & ref["expected"].notna()]
    expected = exp.groupby("species")["expected"].agg(["min", "max"])
    if (expected["max"] - expected["min"]).abs().max() > 1e-12:
        raise SystemExit("fig3: the expected composition differs between the Zymo samples")
    expected = expected["max"]
    truth = list(expected.index)
    spikes = sorted(set(ref.loc[ref["in_truth"] & ref["expected"].isna(), "species"]))
    fp = ref[~ref["in_truth"]]
    bacillus_fp = (fp[fp["species"].str.startswith("Bacillus ")].groupby("species")["estimate"].sum()
                   .sort_values(ascending=False).index.tolist())
    other_fp = sorted(set(fp["species"]) - set(bacillus_fp))
    truth_bacillus = [x for x in truth if x.startswith("Bacillus ")]
    truth_rest = [x for x in truth if x not in truth_bacillus]

    # segments: (key, label, italic, colour, hatch)
    palette = list(style.OKABE_ITO.values())
    palette = [palette[i] for i in (5, 1, 2, 3, 4, 6, 7, 0)]  # blue, orange, sky blue, ...
    truth_colours = dict(zip(truth_bacillus + truth_rest, palette))
    greys = np.linspace(0.30, 0.78, max(len(bacillus_fp), 1))
    segments = [(x, x, True, truth_colours[x], None) for x in truth_bacillus]
    segments += [(x, x, True, str(g), "////" if i % 2 else None)
                 for i, (x, g) in enumerate(zip(bacillus_fp, greys))]
    segments += [(x, x, True, truth_colours[x], None) for x in truth_rest]
    segments += [("__spikes__", "ATCC spike-ins (MSPlus;\npresence only, no expected share)", False,
                  "white", "//////")]
    segments += [("__other_fp__", "other false-positive species", False, "#E3E3E3", None)]

    def composition(df, value):
        v = df.set_index("species")[value]
        out = {}
        for key, *_ in segments:
            if key == "__spikes__":
                out[key] = v.reindex(spikes).fillna(0).sum()
            elif key == "__other_fp__":
                out[key] = v.reindex(other_fp).fillna(0).sum()
            else:
                out[key] = v.get(key, 0.0)
        return out

    bars = [("Expected (Zymo 16S composition)", None, {k: expected.get(k, 0.0) for k, *_ in segments})]
    for sid, label in samples:
        bars.append((label, sid, composition(ref[ref["sample"] == sid], "estimate")))

    fig = style.figure(style.DOUBLE, 92)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.6, 1])
    ax = fig.add_subplot(gs[0, 0])
    lax = fig.add_subplot(gs[0, 1])
    lax.axis("off")
    ypos = []
    y = 0.0
    for i, _ in enumerate(bars):
        ypos.append(y)
        y -= 1.0 + (0.5 if i in (0, 2) else 0.0)
    rows = []
    for (label, sid, comp), yy in zip(bars, ypos):
        left = 0.0
        for key, seg_label, _, colour, hatch in segments:
            w = comp[key]
            if w > 0:
                ax.barh(yy, w, left=left, height=0.78, color=colour, hatch=hatch,
                        edgecolor="#777777" if key == "__spikes__" else "white",
                        hatchcolor="#777777" if key == "__spikes__" else "white",
                        linewidth=0.25)
            left += w
            rows.append({"bar": label, "sample": sid or "expected", "value_type":
                         "expected" if sid is None else "estimate", "segment": seg_label.replace("\n", " "),
                         "species": ";".join(spikes if key == "__spikes__" else other_fp
                                             if key == "__other_fp__" else [key]),
                         "value": w})
    ax.set_yticks(ypos, [b[0] for b in bars])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(ypos[-1] - 0.6, 0.6)
    ax.set_xlabel("Relative abundance")
    ax.spines["left"].set_visible(False)

    def handle(colour, hatch, key):
        return Patch(facecolor=colour, hatch=hatch, edgecolor="#777777" if key == "__spikes__"
                     else "white", linewidth=0.25)

    exp_segments = [sg for sg in segments if sg[0] in truth]
    other_segments = [sg for sg in segments if sg[0] not in truth]
    leg1 = lax.legend([handle(c, h, k) for k, _, _, c, h in exp_segments],
                      [lbl for _, lbl, *_ in exp_segments], title="Expected species",
                      loc="upper left", bbox_to_anchor=(0, 1.0), alignment="left")
    leg2 = lax.legend([handle(c, h, k) for k, _, _, c, h in other_segments],
                      [lbl for _, lbl, *_ in other_segments], title="Spike-ins and false positives",
                      loc="upper left", bbox_to_anchor=(0, 0.50), alignment="left")
    lax.add_artist(leg1)
    for leg, segs in ((leg1, exp_segments), (leg2, other_segments)):
        for t, sg in zip(leg.get_texts(), segs):
            if sg[2]:
                t.set_fontstyle("italic")
    style.panel_labels(fig, [])
    data = pd.DataFrame(rows)
    data["source_file"] = name
    return fig, data


# ─── Figure 4: near-neighbour pairs and detection limit ─────────────────────

def fig4_pairs_lod():
    pairs_name = "pairs_ties.tsv"
    lod_name = "fig_lod.tsv"
    p = read(pairs_name)
    p = p[(p["N"] == 3) & (p["depth"] == "10k")].copy()
    p["realised_share"] = p["reads_a"] / (p["reads_a"] + p["reads_b"])
    p["vb_error"] = p["vb_share_a_in_pair"] - p["realised_share"]
    pair_order = list(dict.fromkeys(p["pair"]))
    combos = [(c, s) for c in ("r1041", "r941") for s in ("even", "log")
              if ((p["chemistry"] == c) & (p["shape"] == s)).any()]
    offs = np.linspace(-0.27, 0.27, len(combos))
    rng = np.random.default_rng(4)

    fig = style.figure(style.DOUBLE, 78)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1])
    axa = fig.add_subplot(gs[0, 0])
    axa.axvline(0, color="#555555", lw=0.6, zorder=0)
    for i, pair in enumerate(pair_order):
        y = len(pair_order) - 1 - i
        for (chem, shape), off in zip(combos, offs):
            d = p[(p["pair"] == pair) & (p["chemistry"] == chem) & (p["shape"] == shape)]
            axa.plot(d["vb_error"], y + off + rng.uniform(-0.035, 0.035, len(d)), ms=3.2, alpha=0.85,
                     **sample_marker(chem, shape, style.NANOVI))
            m = d["vb_error"].mean()
            axa.plot([m, m], [y + off - 0.08, y + off + 0.08], color="black", lw=1.0)
        if i:
            axa.axhline(y + 0.5, color="#E5E5E5", lw=0.4, zorder=0)
    labels = [" / ".join(short_species(x.strip()) for x in pair.split("/")) for pair in pair_order]
    axa.set_yticks(range(len(pair_order) - 1, -1, -1), labels, fontstyle="italic")
    axa.tick_params(axis="y", length=0)
    axa.set_ylim(-0.55, len(pair_order) - 0.45)
    span = p["vb_error"].max() - min(p["vb_error"].min(), 0)
    axa.set_xlim(min(p["vb_error"].min(), 0) - 0.06 * span, p["vb_error"].max() + 0.06 * span)
    axa.xaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: f"{round(v, 6) + 0.0:.1f}".replace("-", "\u2212")))
    axa.set_xlabel("VB share of first member \u2212 realised share")
    handles = [Line2D([], [], ms=3.4, **sample_marker(c, s, style.NANOVI),
                      label=f"{style.SHAPE_LABELS[s]}, {style.CHEMISTRY_LABELS[c]}")
               for c, s in combos]
    handles.append(Line2D([], [], color="black", lw=1.0, label="mean over seeds"))
    handles_a = handles

    lod = read(lod_name)
    arm = "exact_N3_primary_emu_vb_f0.001"
    if arm not in set(lod["arm"]):
        raise SystemExit(f"fig4: arm {arm} not in {lod_name}")
    lod = lod[(lod["arm"] == arm) & (~lod["pair_member"].astype(bool))].copy()
    bins = list(dict.fromkeys(lod["reads_bin"]))
    lod["x"] = lod["reads_bin"].map({b: i for i, b in enumerate(bins)})
    axb = fig.add_subplot(gs[0, 1])
    chem_styles = {"r1041": "-", "r941": "--"}
    chem_dodge = {"r1041": -0.07, "r941": 0.07}
    depths = sorted(lod["depth"].unique())
    for depth in depths:
        for chem, ls in chem_styles.items():
            d = lod[(lod["depth"] == depth) & (lod["chemistry"] == chem)].sort_values("x")
            if d.empty:
                continue
            axb.plot(d["x"] + chem_dodge[chem], d["detected"], color=style.DEPTH_COLOURS[depth], ls=ls, lw=1.0,
                     marker=style.DEPTH_MARKERS[depth], ms=3.0,
                     markerfacecolor=style.DEPTH_COLOURS[depth] if chem == "r1041" else "white")
    axb.set_xticks(range(len(bins)), [b.replace("-", "\u2013") for b in bins], rotation=35,
                   ha="right", rotation_mode="anchor")
    axb.set_xlim(-0.4, len(bins) - 0.6)
    axb.set_ylim(-0.04, 1.06)
    axb.set_xlabel("Reads of the species in the sample")
    axb.set_ylabel("Share of species detected")
    handles = [Line2D([], [], color=style.DEPTH_COLOURS[d], marker=style.DEPTH_MARKERS[d], ms=3.0,
                      lw=1.0, label=f"{depth_label(d)} reads per sample") for d in depths]
    handles += [Line2D([], [], color="#555555", ls=ls, lw=1.0, marker="o", ms=3.0,
                       markerfacecolor="#555555" if c == "r1041" else "white",
                       label=style.CHEMISTRY_LABELS[c]) for c, ls in chem_styles.items()]
    fig.legend(handles=handles_a + handles, loc="outside lower center", ncol=5,
               handletextpad=0.4, columnspacing=1.4)
    style.panel_labels(fig, [(axa, "a"), (axb, "b")])

    da = p[["sample", "shape", "chemistry", "depth", "N", "pair", "a", "b", "reads_a", "reads_b",
            "vb_share_a_in_pair", "realised_share", "vb_error"]].copy()
    da.insert(0, "panel", "a")
    da["source_file"] = pairs_name
    db = lod[["arm", "pair_member", "depth", "chemistry", "reads_bin", "n", "detected"]].copy()
    db.insert(0, "panel", "b")
    db["source_file"] = lod_name
    return fig, pd.concat([da, db], ignore_index=True)


# ─── Figure 5: credible intervals ───────────────────────────────────────────

METHODS = {  # method value -> (label, colour, marker, x offset)
    "mean-field": ("mean-field (shipped)", style.NANOVI, "o", -0.07),
    "gibbs": ("collapsed Gibbs", "#555555", "s", 0.0),
    "bootstrap": ("read bootstrap", "#A0A0A0", "^", 0.07),
}


def fig5_intervals():
    cov_name = "fig_intervals.tsv"
    srf_name = "intervals_coverage_by_shared_read_fraction.tsv"
    auc_name = "intervals_warning_auc.tsv"
    cov = read(cov_name)
    srf = read(srf_name)
    auc = read(auc_name)
    missing = set(cov["method"]) - set(METHODS)
    if missing:
        raise SystemExit(f"fig5: unknown methods {missing} in {cov_name}")
    depths = sorted(cov["depth"].unique())
    xpos = {d: i for i, d in enumerate(depths)}

    fig = style.figure(style.DOUBLE, 66)
    gs = fig.add_gridspec(1, 2)
    axa = fig.add_subplot(gs[0, 0])
    axa.axhline(0.95, color="#555555", lw=0.6, ls=":", zorder=0)
    for method, (label, colour, marker, off) in METHODS.items():
        for pm, ls in ((False, "-"), (True, "--")):
            d = cov[(cov["method"] == method) & (cov["pair_member"].astype(bool) == pm)].sort_values("depth")
            axa.plot(d["depth"].map(xpos) + off, d["coverage_design"], color=colour, ls=ls,
                     marker=marker, ms=3.2, lw=0.9, markerfacecolor=colour if not pm else "white")
    axa.set_xticks(range(len(depths)), [f"{depth_label(d)}" for d in depths])
    axa.set_xlim(-0.35, len(depths) - 0.65)
    axa.set_ylim(0, 1.0)
    axa.set_xlabel("Reads per sample")
    axa.set_ylabel("Coverage of 95% intervals")
    handles = [Line2D([], [], color=c, marker=m, ms=3.2, lw=0.9, label=lbl)
               for lbl, c, m, _ in METHODS.values()]
    handles += [Line2D([], [], color="#555555", ls="-", lw=0.9, label="species outside the pairs"),
                Line2D([], [], color="#555555", ls="--", lw=0.9, marker="o", ms=3.2,
                       markerfacecolor="white", label="pair members"),
                Line2D([], [], color="#555555", ls=":", lw=0.6, label="nominal 0.95")]
    handles_a = handles

    axb = fig.add_subplot(gs[0, 1])
    axb.axhline(0.95, color="#555555", lw=0.6, ls=":", zorder=0)
    mf = srf[srf["method"] == "mean-field"].copy()
    bins = list(dict.fromkeys(mf["srf_bin"]))
    mf["x"] = mf["srf_bin"].map({b: i for i, b in enumerate(bins)})
    auc_mf = auc[auc["method"] == "mean-field"].set_index("depth")["auc_shared_read_fraction"]
    handles = []
    for depth in sorted(mf["depth"].unique()):
        d = mf[mf["depth"] == depth].sort_values("x")
        axb.plot(d["x"], d["coverage_design"], color=style.DEPTH_COLOURS[depth], lw=1.0,
                 marker=style.DEPTH_MARKERS[depth], ms=3.0)
        handles.append(Line2D([], [], color=style.DEPTH_COLOURS[depth], marker=style.DEPTH_MARKERS[depth],
                              ms=3.0, lw=1.0,
                              label=f"{depth_label(depth)} reads, mean-field (AUC {auc_mf[depth]:.2f})"))

    def pretty(interval):
        lo, hi = (float(x) for x in interval.strip("(]").split(","))
        return f"{max(lo, 0):g}–{hi:g}"

    axb.set_xticks(range(len(bins)), [pretty(b) for b in bins])
    axb.set_xlim(-0.35, len(bins) - 0.65)
    axb.set_ylim(0, 1.0)
    axb.set_xlabel("shared_read_fraction of the species")
    axb.set_ylabel("Coverage of 95% intervals")
    fig.legend(handles=handles_a + handles, loc="outside lower center", ncol=3,
               columnspacing=1.6)
    style.panel_labels(fig, [(axa, "a"), (axb, "b")])

    da = cov[["method", "depth", "pair_member", "n", "coverage_design"]].copy()
    da.insert(0, "panel", "a")
    da["source_file"] = cov_name
    db = mf[["method", "depth", "srf_bin", "n", "coverage_design"]].copy()
    db.insert(0, "panel", "b")
    db["source_file"] = srf_name
    dc = auc[auc["method"] == "mean-field"][["method", "depth", "n", "misses",
                                             "auc_shared_read_fraction"]].copy()
    dc.insert(0, "panel", "b")
    dc["source_file"] = auc_name
    return fig, pd.concat([da, db, dc], ignore_index=True)


# ─── Figure 6: clade exclusion ──────────────────────────────────────────────

FAMILIES = {  # family value -> (label, colour)
    "none": ("no filter", style.OKABE_ITO["black"]),
    "pid": ("identity filter", style.OKABE_ITO["vermillion"]),
    "len": ("aligned-length filter", style.OKABE_ITO["bluish_green"]),
    "unk": ("unknown component", style.OKABE_ITO["reddish_purple"]),
}


def arm_threshold(arm):
    m = re.search(r"(\d+(?:\.\d+)?)$", arm)
    return float(m.group(1)) if m else 0.0


def arm_label(family, arm):
    t = arm_threshold(arm)
    if family == "pid":
        return f"identity ≥ {t:g}%"
    if family == "len":
        return f"length ≥ {t:,.0f}"
    if family == "unk":
        return f"δ {t:g}"
    return arm


def fig6_clade():
    name = "fig_clade.tsv"
    c = read(name)
    missing = set(c["family"]) - set(FAMILIES)
    if missing:
        raise SystemExit(f"fig6: unknown families {missing} in {name}")
    chem_markers = {"r1041": "o", "r941": "s"}

    fig = style.figure(style.SINGLE, 80)
    ax = fig.add_subplot(1, 1, 1)
    lin = 1e-3
    for a in ("x", "y"):
        getattr(ax, f"set_{a}scale")("symlog", linthresh=lin, linscale=0.4)
    top = max(c["indb_removed"].max(), c["excluded_removed"].max()) * 1.25
    ax.plot([0, top], [0, top], color="#999999", lw=0.6, ls="--", zorder=0)
    notes = []
    for fam, (label, colour) in FAMILIES.items():
        for chem, marker in chem_markers.items():
            d = c[(c["family"] == fam) & (c["chemistry"] == chem)].copy()
            if d.empty:
                continue
            d["t"] = d["arm"].map(arm_threshold)
            d = d.sort_values("t")
            ax.plot(d["indb_removed"], d["excluded_removed"], color=colour, lw=0.6, marker=marker,
                    ms=3.4, markerfacecolor=colour if chem == "r1041" else "white",
                    markeredgecolor=colour, markeredgewidth=0.7, alpha=0.9,
                    zorder=4 if fam == "none" else 2)
            # label the most aggressive arm of each filter
            if fam in ("pid", "len"):
                ext = d.iloc[-1]
            elif fam == "unk":
                ext = d.iloc[0]
            else:
                continue
            notes.append((ext["excluded_removed"], ext["indb_removed"],
                          f"{arm_label(fam, ext['arm'])}, {style.CHEMISTRY_LABELS[chem]}", colour,
                          fam == "len"))
    # the most aggressive arm of each filter, labelled beside its point
    # (left for points near the top, below-right for the aligned-length arms)
    for yv, xv, text, colour, below in notes:
        ax.annotate(text, (xv, yv), xytext=(4, -5) if below else (-5, 0), textcoords="offset points",
                    ha="left" if below else "right", va="top" if below else "center",
                    fontsize=6, color=colour, zorder=5,
                    bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none",
                              alpha=0.85))
    ticks = [0, 1e-3, 1e-2, 1e-1, 1]
    for a in ("x", "y"):
        getattr(ax, f"set_{a}ticks")(ticks, ["0", "0.001", "0.01", "0.1", "1"])
        getattr(ax, f"set_{a}lim")(-lin * 0.3, top)
    ax.set_xlabel("Share of in-database species' reads removed")
    ax.set_ylabel("Share of excluded species' reads removed")
    handles = [Line2D([], [], color=col, lw=0.6, marker="o", ms=3.4, label=lbl)
               for lbl, col in FAMILIES.values()]
    handles += [Line2D([], [], color="#555555", ls="none", marker=m, ms=3.4,
                       markerfacecolor="#555555" if ch == "r1041" else "white",
                       label=style.CHEMISTRY_LABELS[ch]) for ch, m in chem_markers.items()]
    handles.append(Line2D([], [], color="#999999", ls="--", lw=0.6, label="y = x"))
    ax.legend(handles=handles, loc="lower right")
    data = c.copy()
    data["family_label"] = data["family"].map(lambda f: FAMILIES[f][0])
    data["source_file"] = name
    return fig, data


# ─── Figure 7: runtime ──────────────────────────────────────────────────────

NANOVI_STEPS = [("sam_to_bam", "NanoVI: BAM conversion", 0.62),
                ("nanovi_likelihood", "NanoVI: likelihood", 0.32),
                ("nanovi_vb", "NanoVI: VB", 0.0)]
SAMPLE_LABELS = {"emu_zymo_SRR10391201": "SRR10391201 (R9.4.1)",
                 "zymo_MSPlus_sup_rep1": "MSPlus sup rep 1 (R10.4.1)"}
ALIGN_COLOUR = "#C8C8C8"
NANOVI_DEFAULT = (3, 21)  # the pipeline's -N and -k defaults (nextflow.config: N, kmer_size), used only for a label


def fig7_runtime():
    name = "fig_timing.tsv"
    t = read(name)
    steps = set(t["step"])
    need = {"align", "emu", "emu_shipped"} | {s for s, *_ in NANOVI_STEPS}
    if need - steps:
        raise SystemExit(f"fig7: steps {need - steps} not in {name}")
    matched = t[t["mode"] == "matched"]
    shipped = t[t["mode"] == "emu_shipped"]
    samples = list(dict.fromkeys(matched["sample"]))
    settings = (matched[["N", "k"]].drop_duplicates().astype(int)
                .sort_values(["N", "k"], ascending=[True, False]).itertuples(index=False))
    settings = [(s.N, s.k) for s in settings]
    nanovi_colours = {s: style.mix(style.NANOVI, "white", m) for s, _, m in NANOVI_STEPS}

    # build the bars top to bottom
    bars = []  # dicts: y, label, kind, sample, N, k, steps {step: cpu_s}
    y = 0.0
    headers = []
    for sid in samples:
        headers.append((y, SAMPLE_LABELS.get(sid, sid)))
        y -= 0.9
        for n, k in settings:
            d = matched[(matched["sample"] == sid) & (matched["N"] == n) & (matched["k"] == k)]
            cpu = d.set_index("step")["cpu_s"]
            tag = f"-N {n}, k {k}"
            default = " (NanoVI default)" if (n, k) == NANOVI_DEFAULT else ""
            bars.append(dict(y=y, label=f"NanoVI, {tag}{default}", kind="nanovi", sample=sid, N=n, k=k,
                             steps={s: cpu[s] for s in ["align"] + [x for x, *_ in NANOVI_STEPS]}))
            y -= 0.85
            bars.append(dict(y=y, label=f"Emu, same SAM, {tag}", kind="emu_gtdb", sample=sid, N=n, k=k,
                             steps={s: cpu[s] for s in ("align", "emu")}))
            y -= 1.2
        e = shipped[shipped["sample"] == sid].iloc[0]
        bars.append(dict(y=y, label=f"Emu as shipped, -N {int(e['N'])}, k {int(e['k'])}", kind="emu_shipped",
                         sample=sid, N=int(e["N"]), k=int(e["k"]), steps={"emu_shipped": e["cpu_s"]}))
        y -= 1.7

    fig = style.figure(style.DOUBLE, 112)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1])
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1], sharey=axa)
    h = 0.72
    xmax = max(sum(b["steps"].values()) for b in bars)
    for b in bars:
        yy = b["y"]
        if b["kind"] == "emu_shipped":
            axa.barh(yy, b["steps"]["emu_shipped"], height=h, color=style.EMU_OWN_DB, edgecolor="none")
            axb.text(0.5, yy, "one run (alignment not separated)", va="center", ha="left",
                     fontsize=6, color="#555555")
        else:
            align = b["steps"]["align"]
            post = {s: v for s, v in b["steps"].items() if s != "align"}
            colour = style.NANOVI if b["kind"] == "nanovi" else style.EMU_GTDB
            axa.barh(yy, align, height=h, color=ALIGN_COLOUR, edgecolor="none")
            axa.barh(yy, sum(post.values()), left=align, height=h, color=colour, edgecolor="none")
            left = 0.0
            for s, v in post.items():
                c = nanovi_colours.get(s, style.EMU_GTDB)
                axb.barh(yy, v, left=left, height=h, color=c, edgecolor="white", linewidth=0.25)
                left += v
        total = sum(b["steps"].values())
        axa.text(total + xmax * 0.01, yy, f"{total:,.0f}", va="center", ha="left", fontsize=6)
    ticks = [b["y"] for b in bars] + [hy for hy, _ in headers]
    labels = [b["label"] for b in bars] + [hl for _, hl in headers]
    axa.set_yticks(ticks, labels)
    n_bars = len(bars)
    for i, tl in enumerate(axa.get_yticklabels()):
        if i >= n_bars:
            tl.set_fontweight("bold")
    axa.tick_params(axis="y", length=0)
    axa.set_ylim(min(b["y"] for b in bars) - 0.7, headers[0][0] + 0.6)
    axa.set_xlim(0, xmax * 1.12)
    axa.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    axa.set_xlabel("CPU time (s)")
    plt.setp(axb.get_yticklabels(), visible=False)
    axb.tick_params(axis="y", length=0)
    post_max = max(sum(v for s, v in b["steps"].items() if s != "align")
                   for b in bars if b["kind"] != "emu_shipped")
    axb.set_xlim(0, post_max * 1.08)
    axb.set_xlabel("CPU time after the alignment (s)")
    for ax in (axa, axb):
        ax.grid(axis="x", color="#E5E5E5", lw=0.4)
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
    handles = [Patch(color=ALIGN_COLOUR, label="alignment (minimap2, exact GTDB r226)"),
               Patch(color=style.NANOVI, label="NanoVI after the alignment (a)")]
    handles += [Patch(color=nanovi_colours[s], label=f"{lbl} (b)") for s, lbl, _ in NANOVI_STEPS]
    handles += [Patch(color=style.EMU_GTDB, label="Emu abundance on the same SAM"),
                Patch(color=style.EMU_OWN_DB, label="Emu as shipped (its own database and minimap2)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3)
    style.panel_labels(fig, [(axa, "a"), (axb, "b")])

    rows = []
    for b in bars:
        for s, v in b["steps"].items():
            rows.append({"bar": b["label"], "sample": b["sample"], "kind": b["kind"], "N": b["N"],
                         "k": b["k"], "step": s, "cpu_s": v,
                         "panel": "a" if s == "align" or b["kind"] == "emu_shipped" else "a,b"})
    data = pd.DataFrame(rows)
    data["source_file"] = name
    return fig, data


# ─── Figure 8: clinical cohorts ─────────────────────────────────────────────

COHORTS = [("vaginal_cl", "vaginal (20)", "o"), ("vaginal_bv", "vaginal BV/control (12)", "^"),
           ("crc", "stool CRC/healthy (117)", "s")]
# (arm_a, arm_b, names) of clinical_summary.tsv, top to bottom
CLINICAL_PAIRS = [
    (("nanovi_emudb", "emu_same_sam", "ncbi"), "NanoVI vs Emu,\nsame alignments"),
    (("nanovi_emudb", "emu", "ncbi"), "NanoVI vs Emu,\neach as shipped (Emu db)"),
    (("emu", "published", "ncbi"), "Emu rerun vs\nauthors' Emu tables"),
    (("nanovi_gtdb", "nanovi_emudb", "gtdb"), "NanoVI, GTDB r226\nvs Emu db"),
]
# clinical arms in the fixed tool colours; Emu on NanoVI's alignments is Emu's colour, open
CLINICAL_ARMS = {
    "nanovi_gtdb": ("NanoVI (GTDB r226)", style.NANOVI, True),
    "nanovi_emudb": ("NanoVI (Emu db)", style.TOOLS[("nanovi", "emu_default")][1], True),
    "emu": ("Emu (Emu db)", style.EMU_OWN_DB, True),
    "emu_same_sam": ("Emu on NanoVI's alignments", style.EMU_OWN_DB, False),
    "wgs": ("Kraken2 + Bracken, shotgun", style.TOOLS[("kraken2_bracken", "ncbi_16s_18s")][1], True),
}
BC_FLOOR = 1e-5  # Bray-Curtis values below this are drawn at it (log axis)


def fig8_clinical():
    pairs_name, samples_name, paper_name = ("clinical_pairs.tsv", "clinical_samples.tsv",
                                            "clinical_crc_paper_species.tsv")
    p = read(pairs_name)
    s = read(samples_name)
    paper = read(paper_name)

    fig = style.figure(style.DOUBLE, 128)
    top, bottom = fig.subfigures(2, 1, height_ratios=[1.2, 1])
    axa, axb = top.subplots(1, 2, sharey=True)
    axc, axd = bottom.subplots(1, 2, width_ratios=[1, 1.15])
    rows = []

    # (a, b) per-sample Bray-Curtis for each comparison, cohorts offset within a row
    rng = np.random.default_rng(8)
    grey = "#4D4D4D"
    for i, ((a, b, names), label) in enumerate(CLINICAL_PAIRS):
        y0 = len(CLINICAL_PAIRS) - 1 - i
        for j, (cohort, _, marker) in enumerate(COHORTS):
            d = p[(p["cohort"] == cohort) & (p["arm_a"] == a) & (p["arm_b"] == b) & (p["names"] == names)]
            if d.empty:
                continue
            y = y0 + 0.25 - 0.25 * j
            jit = rng.uniform(-0.07, 0.07, len(d))
            for ax, metric in ((axa, "species_bray_curtis"), (axb, "genus_bray_curtis")):
                x = d[metric].clip(lower=BC_FLOOR)
                ax.scatter(x, y + jit, s=7, marker=marker, facecolors="white", edgecolors=grey,
                           linewidths=0.5, alpha=0.9, zorder=2)
                med = d[metric].median()
                ax.plot([max(med, BC_FLOOR)] * 2, [y - 0.11, y + 0.11], color="black", lw=1.1, zorder=3,
                        path_effects=[pe.Stroke(linewidth=2.6, foreground="white"), pe.Normal()])
            rows.append(with_source(d[["cohort", "run", "sample", "arm_a", "arm_b", "names",
                                       "species_bray_curtis", "genus_bray_curtis"]],
                                    pairs_name, panel="a,b", comparison=label.replace("\n", " ")))
    axa.set_yticks(range(len(CLINICAL_PAIRS) - 1, -1, -1), [l for _, l in CLINICAL_PAIRS])
    axa.set_ylim(-0.55, len(CLINICAL_PAIRS) - 0.45)
    for ax, what in ((axa, "species"), (axb, "genus")):
        ax.set_xscale("log")
        ax.set_xlim(BC_FLOOR / 1.6, 1.3)
        ticks = [BC_FLOOR, 1e-4, 1e-3, 1e-2, 0.1, 1]
        ax.set_xticks(ticks, ["≤10$^{-5}$", "10$^{-4}$", "0.001", "0.01", "0.1", "1"])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlabel(f"Bray–Curtis dissimilarity per sample, {what} level")
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color="#E5E5E5", lw=0.4, zorder=0)
        ax.set_axisbelow(True)
    plt.setp(axb.get_yticklabels(), visible=False)
    handles = [Line2D([], [], marker=m, ls="none", markerfacecolor="white", markeredgecolor=grey, ms=3.6,
                      label=lab) for _, lab, m in COHORTS]
    handles.append(Line2D([], [], color="black", lw=1.1, label="median"))
    top.legend(handles=handles, loc="outside lower center", ncol=4)

    # (c) Lactobacillus share in control and BV swabs
    bv = s[s["cohort"] == "vaginal_bv"]
    arms = ["nanovi_gtdb", "nanovi_emudb", "emu", "wgs"]
    for gi, group in enumerate(("control", "BV")):
        g = bv[bv["group"] == group]
        for ai, arm in enumerate(arms):
            label, colour, filled = CLINICAL_ARMS[arm]
            x = gi * (len(arms) + 1) + ai
            v = g[f"{arm}_lactobacillus"]
            axc.scatter(np.full(len(v), x) + rng.uniform(-0.15, 0.15, len(v)), v, s=9, color=colour,
                        edgecolors=style.edge_for(colour), linewidths=0.4, zorder=2)
            axc.plot([x - 0.3, x + 0.3], [v.median()] * 2, color="black", lw=1.1, zorder=3,
                     path_effects=[pe.Stroke(linewidth=2.6, foreground="white"), pe.Normal()])
            rows.append(with_source(pd.DataFrame({"sample": g["sample"], "group": group, "arm": arm,
                                                  "lactobacillus": v}), samples_name, panel="c"))
    axc.set_xticks([1.5, 1.5 + len(arms) + 1], ["control (6)", "bacterial vaginosis (6)"])
    axc.tick_params(axis="x", length=0)
    axc.set_ylim(-0.03, 1.03)
    axc.set_ylabel("$\\it{Lactobacillus}$ share")

    # (d) the CRC paper's four species: AUC (CRC above healthy) per arm
    species = list(dict.fromkeys(paper["paper_species"]))
    arms_d = ["nanovi_gtdb", "nanovi_emudb", "emu", "emu_same_sam"]
    for si, sp in enumerate(species):
        for ai, arm in enumerate(arms_d):
            r = paper[(paper["paper_species"] == sp) & (paper["arm"] == arm)].iloc[0]
            label, colour, filled = CLINICAL_ARMS[arm]
            y = len(species) - 1 - si + 0.27 - 0.18 * ai
            if r["tested"]:
                axd.plot(r["auc_crc"], y, marker="o", ms=3.8, ls="none", markerfacecolor=colour if filled else "white",
                         markeredgecolor=style.edge_for(colour), markeredgewidth=0.7, zorder=3)
            rows.append(with_source(r.to_frame().T, paper_name, panel="d"))
        r0 = paper[paper["paper_species"] == sp]
        if not r0["tested"].any():
            axd.text(0.5, len(species) - 1 - si, f"in ≤{r0['prevalence'].max():.0%} of samples, not tested",
                     ha="center", va="center", fontsize=6, color=grey)
    axd.axvline(0.5, color="#999999", lw=0.6, ls=":", zorder=1)
    axd.set_yticks(range(len(species) - 1, -1, -1), [short_species(sp) for sp in species])
    style.italic_ticklabels(axd.get_yticklabels())
    axd.set_ylim(-0.6, len(species) - 0.4)
    axd.set_xlim(0.2, 0.8)
    axd.set_xlabel("AUC, CRC above healthy (Mann–Whitney)")
    axd.tick_params(axis="y", length=0)
    bottom.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3.8, markeredgewidth=0.7,
                                  markerfacecolor=CLINICAL_ARMS[a][1] if CLINICAL_ARMS[a][2] else "white",
                                  markeredgecolor=CLINICAL_ARMS[a][1], label=CLINICAL_ARMS[a][0])
                           for a in arms_d + ["wgs"]],
                  loc="outside lower center", ncol=3)
    style.panel_labels(fig, [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")])

    data = pd.concat(rows, ignore_index=True)
    return fig, data


FIGURES = {
    "fig1_database": fig1_database,
    "fig2_comparators": fig2_comparators,
    "fig3_zymo_composition": fig3_zymo_composition,
    "fig4_pairs_lod": fig4_pairs_lod,
    "fig5_intervals": fig5_intervals,
    "fig6_clade": fig6_clade,
    "fig7_runtime": fig7_runtime,
    "fig8_clinical": fig8_clinical,
}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(REPO / "internal_docs" / "manuscript" / "figures"),
                    help="output folder (default: %(default)s)")
    ap.add_argument("--only", nargs="+", metavar="FIG",
                    help="figures to draw, by full name or prefix (fig3, fig3_zymo_composition)")
    args = ap.parse_args(argv)
    names = list(FIGURES)
    if args.only:
        wanted = []
        for o in args.only:
            hits = [n for n in names if n == o or n.startswith(o + "_")]
            if not hits:
                ap.error(f"no figure {o!r}; choose from {', '.join(names)}")
            wanted += hits
        names = list(dict.fromkeys(wanted))
    font = style.apply_style()
    print(f"font: {font}")
    for n in names:
        fig, data = FIGURES[n]()
        for path in style.save(fig, n, data, args.out):
            print(path)


if __name__ == "__main__":
    main()
