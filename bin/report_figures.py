"""Static figures of the run report (--report; F60, D42), drawn with matplotlib in the shared style of
figure_style.py: read fate, composition at species and genus, a species heatmap above HEATMAP_ABOVE samples, the
flagged species, the ambiguity groups and the cell-level table. A figure with a sample axis is split into
consecutive parts of at most PART_MAX samples (D42), each written as <key>[_<part>].{pdf,svg,png,csv}. draw_all
returns one record per figure: its title and caption (the model numbers them) and, for the read fate and
composition figures, the data of their interactive version (report_interactive.py), which keeps every sample.
"""

import numpy as np
import pandas as pd

import figure_style as fs
from report_model import FLAG_SHARED, HEATMAP_ABOVE, HEATMAP_TAXA, MAX_TAXA, TOP_PER_SAMPLE, split_parts

FATE = [("pct_species", "reported species"), ("pct_below_detection", "Below detection"),
        ("pct_unassigned", "Unassigned")]
FATE_COLORS = [fs.CATEGORY[0], fs.CATEGORY[1], "#9A9A9A"]
RAW_COLOR = "#D0D0D0"
PASSED_COLOR = fs.CATEGORY[2]
FLAG_COLORS = {"not flagged": "#BBBBBB", f"shared_read_fraction >= {FLAG_SHARED:g}": fs.LOW_DEPTH,
               "identical 16S": fs.CATEGORY[0], "both": fs.CATEGORY[3]}
LEVEL_COLORS = {"species": fs.CATEGORY[0], "genus": fs.CATEGORY[1], "family": fs.CATEGORY[2]}


def has_menu(groups):
    """The interactive version gets a group menu when the samples have two or more groups."""
    return len({g for g in groups.values() if g}) >= 2


def row_mm(n):
    """Height of one sample's bar: roomy for a few samples, a 6 pt label's height for many."""
    return 4.0 if n <= 14 else 2.6


def part_keys(key, parts):
    return [key] if len(parts) == 1 else [f"{key}_{i + 1}" for i in range(len(parts))]


def italic_labels(texts, plain=("Other",)):
    for t in texts:
        if t.get_text() not in plain and not t.get_text().startswith("(no "):
            t.set_fontstyle("italic")


def sample_axis(ax, samples, groups):
    """Samples top to bottom on the y axis, labelled at 6 pt, with group brackets on the right."""
    ax.set_ylim(len(samples) - 0.5, -0.5)
    fs.sample_ticks(ax, samples, (), "y", fs.SMALL_PT)
    ax.tick_params(axis="y", length=0)
    g = [groups.get(s, "") for s in samples]
    if any(g):
        fs.group_brackets(ax, samples, g, "y", separator=None, offset=0.015)


# ---- read fate

def read_fate(data, outdir):
    import matplotlib.pyplot as plt
    fate = data["fate"].set_index("sample")
    samples, groups = data["samples"], data["group"]
    has_pct = fate["pct_species"].notna().all() if "pct_species" in fate else False
    has_raw = fate["raw_reads"].notna().all() and fate["passed_reads"].notna().all()
    if not has_pct and not has_raw:
        return None
    parts = split_parts(samples)
    keys = part_keys("read_fate", parts)
    out = []
    lo = hi = None
    if has_raw:
        vals = pd.concat([fate["raw_reads"], fate["passed_reads"]]).astype(float)
        vals = vals[vals > 0]
        lo = 10 ** np.floor(np.log10(vals.min() / 1.5)) if len(vals) else 1
        hi = 10 ** np.ceil(np.log10(max(vals.max(), 10))) if len(vals) else 10
    for key, part in zip(keys, parts):
        n = len(part)
        h = 24 + n * row_mm(n)
        ncols = int(has_pct) + int(has_raw)
        fig = plt.figure(figsize=fs.figsize("double" if ncols == 2 else "single", h), layout="constrained")
        axes = fig.subplots(1, ncols, sharey=True, squeeze=False)[0]
        y = np.arange(n)
        rows = []
        i = 0
        if has_pct:
            ax = axes[i]
            i += 1
            left = np.zeros(n)
            for (col, label), color in zip(FATE, FATE_COLORS):
                v = fate.loc[part, col].astype(float).to_numpy()
                ax.barh(y, v, left=left, color=color, height=0.8, label=label, lw=0)
                left += v
                rows += [{"sample": s, "fate": label, "pct_of_passed_reads": float(x)} for s, x in zip(part, v)]
            ax.set_xlim(0, 100)
            ax.set_xlabel("% of the reads that passed fastplong")
            sample_axis(ax, part, groups if ncols == 1 else {})
        if has_raw:
            ax = axes[i]
            raw = fate.loc[part, "raw_reads"].astype(float).to_numpy()
            passed = fate.loc[part, "passed_reads"].astype(float).to_numpy()
            ax.barh(y, raw, color=RAW_COLOR, height=0.8, label="raw reads", lw=0)
            ax.barh(y, passed, color=PASSED_COLOR, height=0.45, label="passed fastplong", lw=0)
            ax.set_xscale("log")
            ax.set_xlim(lo, hi)
            fs.count_axis(ax.xaxis)
            ax.set_xlabel("reads (log scale)")
            if not has_pct:
                sample_axis(ax, part, groups)
            else:
                ax.tick_params(axis="y", length=0)
                g = [groups.get(s, "") for s in part]
                if any(g):
                    fs.group_brackets(ax, part, g, "y", separator=None, offset=0.015)
            rows += [{"sample": s, "fate": "raw reads", "reads": float(r)} for s, r in zip(part, raw)]
            rows += [{"sample": s, "fate": "passed fastplong", "reads": float(p)} for s, p in zip(part, passed)]
        handles = [hd for a in axes for hd in a.get_legend_handles_labels()[0]]
        fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False, handlelength=1,
                   columnspacing=1)
        src = pd.DataFrame(rows)
        if any(groups.values()):
            src.insert(1, "group", [groups.get(s, "") for s in src["sample"]])
        fs.save(fig, outdir, key, src)
        out.append({"key": key, "width_mm": fs.WIDTH_MM["double" if ncols == 2 else "single"],
                    "height_mm": round(min(h, fs.MAX_H_MM), 1)})
    inter = None
    if has_pct or has_raw:
        cols = ["sample"] + [c for c, _ in FATE if has_pct] + (["raw_reads", "passed_reads"] if has_raw else [])
        frame = data["fate"][cols].copy()
        frame["group"] = [groups.get(s, "") for s in frame["sample"]]
        inter = {"kind": "read_fate", "frame": frame, "samples": samples, "has_pct": bool(has_pct),
                 "has_raw": bool(has_raw), "domain": [lo, hi] if has_raw else None, "menu": has_menu(groups)}
    caption = ("Left: where the reads that passed fastplong's length filter went: reported species, Below detection "
               "(species pruned by the detection rules) and Unassigned (no alignment). " if has_pct else "") + \
              ("Right: raw reads and reads that passed fastplong (log scale)." if has_raw else "")
    if len(parts) > 1:
        caption += f" The {len(samples)} samples are split into {len(parts)} parts in sample order."
    return {"key": "read_fate", "section": "overview", "title": "Read fate per sample", "caption": caption.strip(),
            "parts": out, "interactive": inter}


# ---- composition

def composition(data, outdir, rank):
    import matplotlib.pyplot as plt
    long, shown = data[f"composition_{rank}"], data[f"shown_{rank}"]
    samples, groups = data["samples"], data["group"]
    taxa = list(shown) + ["Other"]
    colors = fs.category_colors(taxa)
    wide = long.pivot(index="sample", columns="taxon", values="abundance").reindex(index=samples, columns=taxa)
    parts = split_parts(samples)
    keys = part_keys(f"composition_{rank}", parts)
    out = []
    for key, part in zip(keys, parts):
        n = len(part)
        h = max(20 + n * row_mm(n), 14 + len(taxa) * 3.0)
        fig = plt.figure(figsize=fs.figsize("double", h), layout="constrained")
        ax = fig.subplots()
        y = np.arange(n)
        left = np.zeros(n)
        for t in taxa:
            v = 100 * wide.loc[part, t].astype(float).to_numpy()
            ax.barh(y, v, left=left, color=colors[t], height=0.8, label=fs.clip_label(t, 45), lw=0)
            left += v
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of the reads on reported species (read-level abundance)")
        sample_axis(ax, part, groups)
        leg = fig.legend(loc="outside right upper", frameon=False, handlelength=1, fontsize=fs.SMALL_PT,
                         title="genus" if rank == "genus" else "species", alignment="left")
        italic_labels(leg.get_texts())
        src = long[long["sample"].isin(part)].copy()
        if any(groups.values()):
            src.insert(1, "group", [groups.get(s, "") for s in src["sample"]])
        fs.save(fig, outdir, key, src)
        out.append({"key": key, "width_mm": fs.WIDTH_MM["double"], "height_mm": round(min(h, fs.MAX_H_MM), 1)})
    frame = long.copy()
    frame["group"] = [groups.get(s, "") for s in frame["sample"]]
    inter = {"kind": "bars", "frame": frame, "samples": samples, "taxa": taxa,
             "colors": [colors[t] for t in taxa], "rank": rank, "menu": has_menu(groups)}
    caption = (f"Read-level relative abundance of the {len(shown)} {'genera' if rank == 'genus' else 'species'} "
               f"among each sample's top {TOP_PER_SAMPLE} (at most {MAX_TAXA}, by mean abundance); Other is the "
               "sum of the rest. Below detection and Unassigned reads are not part of these shares.")
    if len(parts) > 1:
        caption += f" The {len(samples)} samples are split into {len(parts)} parts in sample order."
    return {"key": f"composition_{rank}", "section": "composition",
            "title": f"Composition per sample ({rank})", "caption": caption, "parts": out, "interactive": inter}


# ---- heatmaps

def heatmap(wide, samples, groups, outdir, key, label_len=45, italic=True, legend="log10 relative abundance"):
    """rows x samples heatmap of log10 abundance, 'not detected' grey, parts of at most PART_MAX samples."""
    import matplotlib.pyplot as plt
    parts = split_parts(samples)
    keys = part_keys(key, parts)
    vals = wide.to_numpy(dtype=float)
    pos = vals[vals > 0]
    vmin = np.floor(np.log10(pos.min())) if pos.size else -3
    vmax = 0.0
    out = []
    longest = max([len(str(s)) for s in samples] + [1])
    for k, part in zip(keys, parts):
        m = wide[part].to_numpy(dtype=float)
        with np.errstate(divide="ignore"):
            logv = np.where(m > 0, np.log10(np.where(m > 0, m, 1)), np.nan)
        h = 30 + wide.shape[0] * 2.6 + min(longest, 30) * 1.1
        fig = plt.figure(figsize=fs.figsize("double", h), layout="constrained")
        ax = fig.subplots()
        im = ax.imshow(np.ma.masked_invalid(logv), aspect="auto", cmap=fs.viridis(), vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_yticks(range(wide.shape[0]), [fs.clip_label(t, label_len) for t in wide.index], fontsize=fs.SMALL_PT)
        if italic:
            italic_labels(ax.get_yticklabels())
        fs.sample_ticks(ax, part, (), "x", fs.sample_layout(len(part))["fs"])
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        g = [groups.get(s, "") for s in part]
        if any(g):
            fs.group_brackets(ax, part, g, "x")
        cb = fig.colorbar(im, ax=ax, location="bottom", shrink=0.35, aspect=30, pad=0.02)
        cb.set_label(legend, fontsize=fs.SMALL_PT)
        cb.ax.tick_params(labelsize=fs.SMALL_PT)
        src = wide[part].reset_index().melt(id_vars=wide.index.name or "index", var_name="sample",
                                            value_name="abundance")
        src = src.rename(columns={wide.index.name or "index": "row"})
        fs.save(fig, outdir, k, src[["sample", "row", "abundance"]])
        out.append({"key": k, "width_mm": fs.WIDTH_MM["double"], "height_mm": round(min(h, fs.MAX_H_MM), 1)})
    return out, parts, [float(vmin), vmax]


def species_heatmap(data, outdir):
    wide = data["heatmap_species"].copy()
    wide.index.name = "species"
    if wide.empty:
        return None
    samples, groups = data["samples"], data["group"]
    out, parts, domain = heatmap(wide, samples, groups, outdir, "heatmap_species")
    frame = wide.reset_index().melt(id_vars="species", var_name="sample", value_name="abundance")
    frame["group"] = [groups.get(s, "") for s in frame["sample"]]
    caption = (f"The {len(wide)} species with the highest mean read-level abundance (at most {HEATMAP_TAXA}; drawn "
               f"above {HEATMAP_ABOVE} samples), log10 scale; grey: not detected.")
    if len(parts) > 1:
        caption += f" The {len(samples)} samples are split into {len(parts)} parts in sample order."
    return {"key": "heatmap_species", "section": "composition", "title": "Species heatmap", "caption": caption,
            "parts": out, "interactive": {"kind": "heatmap", "frame": frame, "samples": samples,
                                          "rows": list(wide.index), "domain": domain, "feature": "species",
                                          "menu": has_menu(groups)}}


def groups_heatmap(data, outdir):
    wide = data.get("groups_wide")
    if wide is None or wide.empty:
        return None
    wide = wide.copy()
    wide.index.name = "members"
    samples = data["samples"]
    out, parts, _ = heatmap(wide, samples, data["group"], outdir, "groups_heatmap", label_len=60)
    caption = (f"Read-level abundance of the ambiguity groups of two or more species (the {len(wide)} member sets "
               f"with the highest mean, at most {HEATMAP_TAXA}), log10 scale; grey: no such group in the sample.")
    if len(parts) > 1:
        caption += f" The {len(samples)} samples are split into {len(parts)} parts in sample order."
    return {"key": "groups_heatmap", "section": "groups", "title": "Ambiguity groups per sample", "caption": caption,
            "parts": out}


# ---- scatter plots

def flagged_species(data, outdir):
    import matplotlib.pyplot as plt
    long = data["flagged_long"]
    d = long[long["shared_read_fraction"].notna() & (long["abundance"] > 0)].copy()
    if d.empty:
        return None
    names = list(FLAG_COLORS)
    d["flag"] = np.select([d["flag_shared"] & d["flag_identical"], d["flag_shared"], d["flag_identical"]],
                          [names[3], names[1], names[2]], default=names[0])
    fig = plt.figure(figsize=fs.figsize("single", 78), layout="constrained")
    ax = fig.subplots()
    for name in names:
        sub = d[d["flag"] == name]
        if len(sub):
            ax.scatter(sub["shared_read_fraction"], sub["abundance"], s=6 if name == names[0] else 12,
                       color=FLAG_COLORS[name], lw=0, label=f"{name} ({len(sub)})", zorder=2 if name != names[0] else 1)
    ax.axvline(FLAG_SHARED, color="#555555", lw=0.5, ls=(0, (3, 2)))
    ax.set_yscale("log")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("shared_read_fraction")
    ax.set_ylabel("read-level abundance (log scale)")
    fig.legend(loc="outside lower center", ncol=2, frameon=False, handlelength=1, markerscale=1.5)
    src = d[["sample", "species", "abundance", "shared_read_fraction", "identical_16s_species", "flag"]]
    fs.save(fig, outdir, "flagged_species", src)
    return {"key": "flagged_species", "section": "flagged", "title": "Flagged species",
            "caption": f"Every reported species in every sample ({len(d)} points): abundance against "
                       f"shared_read_fraction, coloured by flag (dashed line: {FLAG_SHARED:g}).",
            "parts": [{"key": "flagged_species", "width_mm": fs.WIDTH_MM["single"], "height_mm": 78}]}


def cell_vs_read(data, outdir):
    import matplotlib.pyplot as plt
    cl = data.get("cell_long")
    if cl is None or cl.empty or "read_abundance" not in cl.columns:
        return None
    d = cl[(cl["read_abundance"] > 0) & (cl["abundance"] > 0)].copy()
    if d.empty:
        return None
    fig = plt.figure(figsize=fs.figsize("single", 80), layout="constrained")
    ax = fig.subplots()
    lo = 10 ** np.floor(np.log10(min(d["read_abundance"].min(), d["abundance"].min())))
    ax.plot([lo, 1], [lo, 1], color="#888888", lw=0.5, ls=(0, (3, 2)), zorder=0)
    for level, color in LEVEL_COLORS.items():
        sub = d[d["copy_number_level"] == level]
        if len(sub):
            ax.scatter(sub["read_abundance"], sub["abundance"], s=12, color=color, lw=0,
                       label=f"copy number of the {level} ({len(sub)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, 1.2)
    ax.set_ylim(lo, 1.2)
    ax.set_xlabel("read-level abundance")
    ax.set_ylabel("cell-level abundance")
    fig.legend(loc="outside lower center", ncol=1, frameon=False, handlelength=1, markerscale=1.5)
    src = d[["sample", "species", "read_abundance", "abundance", "copy_number", "copy_number_level"]]
    fs.save(fig, outdir, "cell_vs_read", src)
    return {"key": "cell_vs_read", "section": "cell", "title": "Cell-level against read-level abundance",
            "caption": f"Each species in each sample ({len(d)} points), log scales; dashed: equal abundances. A "
                       "species above the line has fewer 16S copies than the sample's average.",
            "parts": [{"key": "cell_vs_read", "width_mm": fs.WIDTH_MM["single"], "height_mm": 80}]}


def draw_all(data, outdir):
    """Every figure the data allow; returns their records."""
    fs.setup()
    records = []
    for r in (read_fate(data, outdir), composition(data, outdir, "species"), composition(data, outdir, "genus"),
              species_heatmap(data, outdir) if "heatmap_species" in data else None,
              flagged_species(data, outdir), groups_heatmap(data, outdir), cell_vs_read(data, outdir)):
        if r is not None:
            records.append(r)
    return records
