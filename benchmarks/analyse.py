#!/usr/bin/env python3
"""Tables and figures from the Phase 7a results (items 3, 5 and 9).

  ablation       one-factor and cell grids: mean and SD over seeds, paired differences
                 against the reference setting
  contamination  c995 and reps against exact: reads whose primary hit is in another
                 genus, and how many of those hit a sequence that GTDB has in more than
                 one genus (cross-genus identical, likely contamination; db_audit.py)
  pairs          near-neighbour pairs: reads whose best scores tie between the pair,
                 and how each estimator splits the pair
  lod            detection rate by realised reads, depth and chemistry
  clade          clade exclusion: where an excluded species' reads go, and how much more
                 often each filter or the unknown component removes them than reads of
                 species in the database
  real           real mocks: each arm by sample group, and paired against the reference
  intervals      credible interval coverage by method, depth, pair membership, abundance
                 and shared_read_fraction, and how well that fraction flags the misses

Every figure (docs/benchmarks/figures/*.png) has its source data in
docs/benchmarks/data/ under the same name.
"""
import argparse
import hashlib
import json
from collections import Counter
from functools import lru_cache

import numpy as np
import pandas as pd
import pysam

from common import BENCH, DOCS_DATA, DOCS_FIG, write_tsv
from evaluate import load_estimate, taxonomy
from likelihood import score_reads
from simulate import PAIRS

RESULTS = BENCH / "results"
METRICS = ["species_l1", "genus_l1", "group_l1", "species_recall", "species_precision",
           "species_fp_count", "species_fp_abundance"]
REFERENCE_ARM = "exact_N3_primary_emu_vb_f0.001"


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    return plt


def save_figure(fig, name, data):
    DOCS_FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(DOCS_FIG / f"{name}.png", bbox_inches="tight")
    write_tsv(data, DOCS_DATA / f"{name}.tsv")


# ─── Ablation ───────────────────────────────────────────────────────────────

def ablation():
    ofat = pd.read_csv(RESULTS / "grid_ofat.tsv", sep="\t")
    keys = ["depth", "chemistry", "shape", "arm"]
    summary = ofat.groupby(keys)[METRICS].agg(["mean", "std"])
    summary.columns = [f"{m}_{s}" for m, s in summary.columns]
    summary["n_samples"] = ofat.groupby(keys).size()
    write_tsv(summary.reset_index(), DOCS_DATA / "ablation_ofat_summary.tsv")

    ref = ofat[ofat["arm"] == REFERENCE_ARM].set_index("sample")
    rows = []
    for arm, t in ofat[ofat["arm"] != REFERENCE_ARM].groupby("arm"):
        t = t.set_index("sample")
        for m in METRICS:
            d = t[m] - ref.loc[t.index, m]
            for (depth, chem), dd in d.groupby([t["depth"], t["chemistry"]]):
                rows.append({"arm": arm, "metric": m, "depth": depth, "chemistry": chem, "n": len(dd),
                             "mean_diff": dd.mean(), "better": int((dd < 0).sum() if "l1" in m or "fp" in m
                                                                   else (dd > 0).sum()),
                             "worse": int((dd > 0).sum() if "l1" in m or "fp" in m else (dd < 0).sum())})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "ablation_ofat_paired.tsv")

    cell = pd.read_csv(RESULTS / "grid_cell.tsv", sep="\t")
    pruned = cell[cell["prune"]]
    for name, keys in (("db_N", ["db", "N"]), ("score_profile", ["score", "profile"]),
                       ("estimator_floor", ["estimator", "floor", "prune"])):
        src = cell if name == "estimator_floor" else pruned
        t = src.groupby(keys)[METRICS].mean().reset_index()
        write_tsv(t, DOCS_DATA / f"ablation_cell_{name}.tsv")

    # Figure: species and genus L1 by database and -N, reference settings otherwise
    plt = _plt()
    base = cell[(cell["profile"] == "primary") & (cell["score"] == "emu") & (cell["estimator"] == "vb")
                & (cell["floor"] == 1e-3) & cell["prune"]]
    data = base[["sample", "db", "N", "species_l1", "genus_l1", "group_l1", "species_fp_abundance"]]
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), sharey=False)
    for ax, m in zip(axes, ["species_l1", "group_l1", "genus_l1"]):
        for i, db in enumerate(["exact", "c995", "reps"]):
            d = data[data["db"] == db].groupby("N")[m]
            ax.errorbar(d.mean().index + (i - 1) * 0.6, d.mean(), yerr=d.std(), marker="o", ms=3, capsize=2,
                        label=db)
        ax.set_xscale("log")
        ax.set_xticks([3, 10, 50], ["3", "10", "50"])
        ax.set_xlabel("minimap2 -N")
        ax.set_title(m.replace("_", " "))
    axes[0].set_ylabel("L1 distance to truth (0-2)")
    axes[0].legend(frameon=False)
    save_figure(fig, "fig_db_N", data)


# ─── Contamination in the database variants ─────────────────────────────────

@lru_cache(maxsize=None)
def cross_genus_md5():
    seqs = pd.read_csv(BENCH / "audit" / "ssu_all" / "sequences.tsv.gz", sep="\t",
                       usecols=["md5", "genus"], dtype=str)
    n = seqs.groupby("md5")["genus"].nunique()
    return frozenset(n[n > 1].index)


@lru_cache(maxsize=None)
def ref_md5(db):
    out = {}
    with pysam.FastxFile(str(BENCH / "db" / db / "species_taxid.fasta")) as f:
        for rec in f:
            out[rec.name] = hashlib.md5(rec.sequence.upper().encode()).hexdigest()
    return out


def contamination(samples, dbs=("exact", "c995", "reps"), N=3):
    cross = cross_genus_md5()
    rows, fp_rows = [], []
    for sample in samples:
        truth = pd.read_csv(BENCH / "sims" / sample / "read_truth.tsv.gz", sep="\t").set_index("read")["species"]
        truth_species = set(truth)
        for db in dbs:
            tax = taxonomy(str(BENCH / "db" / db))
            md5 = ref_md5(db)
            hits = []
            with pysam.AlignmentFile(str(BENCH / "runs" / sample / "aln" / f"{db}_N{N}.bam"), "rb") as bam:
                for a in bam:
                    if a.is_unmapped or a.is_secondary or a.is_supplementary:
                        continue
                    sp = tax.loc[int(a.reference_name.split(":")[0]), "species"]
                    true = truth.get(a.query_name)
                    hits.append((true, sp, md5[a.reference_name] in cross))
            h = pd.DataFrame(hits, columns=["true", "hit", "cross_genus_seq"])
            h["wrong_genus"] = h["true"].str.split().str[0] != h["hit"].str.split().str[0]
            rows.append({"sample": sample, "db": db, "primary_reads": len(h),
                         "wrong_genus_reads": int(h["wrong_genus"].sum()),
                         "wrong_genus_cross_genus_seq": int((h["wrong_genus"] & h["cross_genus_seq"]).sum()),
                         "right_genus_cross_genus_seq": int((~h["wrong_genus"] & h["cross_genus_seq"]).sum())})
            est, _ = load_estimate(BENCH / "runs" / sample / "est" / f"{db}_N{N}_primary_emu_vb_f0.001.json",
                                   BENCH / "db" / db)
            by_hit = h.groupby("hit")["cross_genus_seq"].mean()
            for sp, ab in est.items():
                if sp in truth_species:
                    continue
                fp_rows.append({"sample": sample, "db": db, "species": sp, "abundance": ab,
                                "primary_reads": int((h["hit"] == sp).sum()),
                                "share_on_cross_genus_seq": float(by_hit.get(sp, np.nan)),
                                "same_genus_as_a_truth_species": sp.split()[0] in {s.split()[0] for s in truth_species}})
    t = pd.DataFrame(rows)
    fp = pd.DataFrame(fp_rows)
    write_tsv(t, DOCS_DATA / "contamination_db_variants_reads.tsv")
    write_tsv(fp.sort_values(["db", "sample", "abundance"], ascending=[True, True, False]),
              DOCS_DATA / "contamination_db_variants_fp.tsv")
    fp["explained"] = fp["share_on_cross_genus_seq"] >= 0.5
    s = fp.groupby(["db", "sample"]).apply(lambda d: pd.Series({
        "fp_species": len(d), "fp_abundance": d["abundance"].sum(),
        "fp_abundance_cross_genus_seq": d.loc[d["explained"], "abundance"].sum(),
        "fp_abundance_same_genus": d.loc[d["same_genus_as_a_truth_species"], "abundance"].sum()}),
        include_groups=False).reset_index()
    write_tsv(s, DOCS_DATA / "contamination_db_variants_summary.tsv")
    return t, s


# ─── Near-neighbour pairs ───────────────────────────────────────────────────

PAIR_ARMS = {"vb": REFERENCE_ARM, "em": "exact_N3_primary_emu_em_f0.001",
             "vb_noprune": "exact_N3_primary_emu_vb_f0_noprune", "per-base": "exact_N3_primary_per-base_vb_f0.001"}
PAIR_CLASSES = ["tie", "favours_true", "favours_other", "true_not_candidate", "other_not_candidate", "neither",
                "unaligned"]


def classify_pair_reads(truth, log_p_rgs, species_of, a, b):
    """Counts of the pair's reads by how their best scores for a and b compare.

    true_not_candidate: the read aligns to the other member but none of its kept
    alignments (-N) is to its own species, so the estimator cannot give it to it.
    """
    c = Counter()
    for read, true in truth.items():
        if true not in (a, b):
            continue
        c[true] += 1
        cand = log_p_rgs.get(read)
        if cand is None:
            c["unaligned"] += 1
            continue
        best = {}
        for tid, score in zip(*cand):
            sp = species_of[tid]
            best[sp] = max(best.get(sp, -np.inf), score)
        other = b if true == a else a
        st, so = best.get(true), best.get(other)
        if st is None and so is None:
            c["neither"] += 1
        elif st is None:
            c["true_not_candidate"] += 1
        elif so is None:
            c["other_not_candidate"] += 1
        elif st == so:
            c["tie"] += 1
        else:
            c["favours_true" if st > so else "favours_other"] += 1
        if max(best.get(a, -np.inf), best.get(b, -np.inf)) < max(best.values()):
            c["outside_pair_best"] += 1
    return c


def pairs(samples):
    """Every one-factor sample on the reference alignments, and the cell samples at -N 3, 10, 50."""
    tax = taxonomy(str(BENCH / "db" / "exact"))
    species_of = tax["species"].to_dict()
    ofat = pd.read_csv(RESULTS / "grid_ofat_species.tsv.gz", sep="\t")
    ofat = ofat[ofat["arm"].isin(PAIR_ARMS.values())].set_index(["sample", "arm", "species"])["estimate"]
    cell = pd.read_csv(RESULTS / "grid_cell_species.tsv.gz", sep="\t")
    cell_samples = set(cell["sample"])
    cell = cell.set_index(["sample", "arm", "species"])["estimate"]
    rows = []
    for sample in samples:
        truth = pd.read_csv(BENCH / "sims" / sample / "read_truth.tsv.gz", sep="\t").set_index("read")["species"]
        for N in (3, 10, 50) if sample in cell_samples else (3,):
            log_p_rgs, _, _ = score_reads(str(BENCH / "runs" / sample / "aln" / f"exact_N{N}.bam"))
            if N == 3:
                arms_est = [(name, arm, ofat) for name, arm in PAIR_ARMS.items()]
            else:
                arms_est = [("vb", f"exact_N{N}_primary_emu_vb_f0.001", cell),
                            ("em", f"exact_N{N}_primary_emu_em_f0.001", cell),
                            ("vb_noprune", f"exact_N{N}_primary_emu_vb_f0_noprune", cell)]
            for a, b in PAIRS:
                c = classify_pair_reads(truth, log_p_rgs, species_of, a, b)
                n = c[a] + c[b]
                row = {"sample": sample, "N": N, "pair": f"{a} / {b}", "a": a, "b": b, "reads": n,
                       "reads_a": c[a], "reads_b": c[b],
                       **{k: c[k] for k in PAIR_CLASSES + ["outside_pair_best"]}}
                for name, arm, est in arms_est:
                    ea, eb = est.get((sample, arm, a), 0.0), est.get((sample, arm, b), 0.0)
                    row[f"{name}_a"], row[f"{name}_b"] = ea, eb
                    row[f"{name}_share_a_in_pair"] = ea / (ea + eb) if ea + eb else np.nan
                rows.append(row)
    t = pd.DataFrame(rows)
    t["realised_share_a_in_pair"] = t["reads_a"] / t["reads"]
    for col, src in (("depth", 2), ("chemistry", 1), ("shape", 0)):
        t.insert(1, col, t["sample"].str.split("_").str[src])
    write_tsv(t, DOCS_DATA / "pairs_ties.tsv")
    return t


# ─── Detection limit ────────────────────────────────────────────────────────

BINS = [0, 2, 5, 10, 20, 50, 100, np.inf]
BIN_LABELS = ["1-2", "3-5", "6-10", "11-20", "21-50", "51-100", ">100"]


def lod():
    """Detection by realised reads. Pair members are reported apart: tied reads from the
    partner can give them an estimate without reads of their own."""
    sp = pd.read_csv(RESULTS / "grid_ofat_species.tsv.gz", sep="\t")
    sp = sp[sp["in_truth"] & sp["arm"].isin([REFERENCE_ARM, "exact_N3_primary_emu_vb_f0_noprune",
                                             "exact_N3_primary_per-base_vb_f0.001"])]
    sp = sp[sp["reads"] > 0].copy()
    paired = {s for pair in PAIRS for s in pair}
    sp["pair_member"] = sp["species"].isin(paired)
    sp["detected"] = sp["estimate"] > 0
    sp["reads_bin"] = pd.cut(sp["reads"], BINS, labels=BIN_LABELS)
    keys = ["arm", "pair_member", "depth", "chemistry"]
    t = (sp.groupby(keys + ["reads_bin"], observed=True)
         .agg(n=("detected", "size"), detected=("detected", "mean"), median_share=("realised", "median"))
         .reset_index())
    write_tsv(t, DOCS_DATA / "lod_detection.tsv")
    # lowest read bin from which every bin up is detected in at least 90% of cases
    rows = []
    for k, d in t.groupby(keys):
        ok = d["detected"].to_numpy() >= 0.9
        first = next((i for i in range(len(d)) if ok[i:].all()), None)
        rows.append({**dict(zip(keys, k)), "lowest_bin_90pct": None if first is None else d["reads_bin"].iloc[first],
                     "median_share_in_bin": None if first is None else d["median_share"].iloc[first]})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "lod_summary.tsv")

    plt = _plt()
    data = t[~t["pair_member"]]
    depths = sorted(data["depth"].unique())
    fig, axes = plt.subplots(1, len(depths), figsize=(3.5 * len(depths), 2.6), sharey=True)
    for ax, depth in zip(axes, depths):
        for (arm, chem), d in data[data["depth"] == depth].groupby(["arm", "chemistry"]):
            label = {"exact_N3_primary_emu_vb_f0.001": "default", "exact_N3_primary_emu_vb_f0_noprune": "no pruning",
                     "exact_N3_primary_per-base_vb_f0.001": "per-base score"}[arm]
            ax.plot(d["reads_bin"].astype(str), d["detected"], marker="o", ms=3,
                    ls="-" if chem == "r1041" else "--", label=f"{label}, {chem}")
        ax.set_title(f"{depth:,} reads")
        ax.set_xlabel("reads of the species")
    axes[0].set_ylabel("fraction detected")
    axes[-1].legend(frameon=False, fontsize=6)
    save_figure(fig, "fig_lod", data)
    return t


# ─── Clade exclusion ────────────────────────────────────────────────────────

def clade():
    """Removal rates: an excluded species' reads that end up unassigned (filters) or in
    the unknown component, against the same for the reads of the other 19 species. A
    useful novelty signal removes the first far more often than the second."""
    from clade_exclusion import db_name
    d = pd.read_csv(RESULTS / "clade_exclusion.tsv", sep="\t")
    d["chemistry"] = d["sample"].str.split("_").str[1]
    d["sibling_present"] = d["sibling_in_sample"].notna()
    d["family"] = d["arm"].str.extract(r"^(none|pid|len|unk)")[0]
    indb = []
    for r in d.itertuples():
        npz = r.arm if not r.arm.startswith("unk") else "unk" + r.arm[3:]
        c = json.loads((BENCH / "runs" / r.sample / "clade" / db_name(r.excluded) / f"{npz}.counts.json").read_text())
        total = c["assigned_count"] + c["unassigned_count"]
        others = total - r.excluded_reads
        unassigned = c["unassigned_count"] - r.to_unassigned * r.excluded_reads
        unknown = r.unknown_abundance * c["assigned_count"] - r.to_unknown * r.excluded_reads
        indb.append((unassigned + unknown) / others)
    d["excluded_removed"] = d["to_unassigned"] + d["to_unknown"]
    d["indb_removed"] = indb
    cols = ["excluded_removed", "indb_removed", "unknown_abundance", "to_unknown", "to_sibling",
            "to_other_community", "to_absent_species", "to_below_detection", "to_unassigned",
            "indb_l1", "indb_fp_abundance", "indb_recall", "indb_precision"]
    t = d.groupby(["chemistry", "sibling_present", "family", "arm"])[cols].mean().reset_index()
    t["removal_ratio"] = t["excluded_removed"] / t["indb_removed"].where(t["indb_removed"] > 0)
    t["n"] = d.groupby(["chemistry", "sibling_present", "family", "arm"]).size().to_numpy()
    write_tsv(t, DOCS_DATA / "clade_exclusion_summary.tsv")

    plt = _plt()
    data = d.groupby(["chemistry", "family", "arm"])[["excluded_removed", "indb_removed"]].mean().reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.2), sharey=True)
    for ax, chem in zip(axes, ["r941", "r1041"]):
        for fam, label in (("pid", "identity floor"), ("len", "aligned-length floor"), ("unk", "unknown component")):
            x = data[(data["chemistry"] == chem) & (data["family"] == fam)].sort_values("indb_removed")
            ax.plot(x["indb_removed"], x["excluded_removed"], marker="o", ms=3, label=label)
        ax.plot([0, 1], [0, 1], color="grey", lw=0.6, ls=":")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_title(chem)
        ax.set_xlabel("in-database reads removed")
    axes[0].set_ylabel("excluded species' reads removed")
    axes[0].legend(frameon=False, fontsize=6)
    save_figure(fig, "fig_clade", data)
    return t


# ─── Real mocks ─────────────────────────────────────────────────────────────

def real():
    d = pd.read_csv(RESULTS / "real_nanovi.tsv", sep="\t")
    d["group"] = d["sample"].str.replace(r"_rep\d+$", "", regex=True)
    cols = [c for c in d.columns if c.startswith(("species_", "genus_", "group_"))] + [
        "assigned_count", "unassigned_count", "below_detection"]
    t = d.groupby(["group", "arm"])[cols].agg(["mean", "std"])
    t.columns = [f"{m}_{s}" for m, s in t.columns]
    t["n_samples"] = d.groupby(["group", "arm"]).size()
    write_tsv(t.reset_index(), DOCS_DATA / "real_arms_summary.tsv")
    ref = d[d["arm"] == "reference"].set_index("sample")
    rows = []
    for arm, x in d[d["arm"] != "reference"].groupby("arm"):
        x = x.set_index("sample")
        for m in ("species_l1_known", "genus_l1_known", "species_fp_abundance", "species_fp_count",
                  "species_recall", "species_precision"):
            diff = x[m] - ref.loc[x.index, m]
            for group, dd in diff.groupby(x["group"]):
                rows.append({"arm": arm, "group": group, "metric": m, "n": len(dd), "mean_diff": dd.mean(),
                             "min_diff": dd.min(), "max_diff": dd.max()})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "real_arms_paired.tsv")

    # per species: estimate against the expected 16S share (reference, N50, emudb)
    from evaluate import load_real_truth
    from names import load as load_crosswalk, to_gtdb
    from real import EMU_DB
    table = load_crosswalk()
    rows = []
    for sample in sorted(d["sample"].unique()):
        truth = load_real_truth(sample)
        for arm, db in (("reference", BENCH / "db" / "exact"), ("N50", BENCH / "db" / "exact"), ("emudb", EMU_DB)):
            path = BENCH / "real" / "runs" / sample / "est" / f"{arm}.json"
            if not path.exists():
                continue
            est = {}
            for k, v in load_estimate(path, db)[0].items():
                k = to_gtdb(k, table) if arm == "emudb" else k
                est[k] = est.get(k, 0.0) + v
            for sp in sorted(set(est) | set(truth)):
                rows.append({"sample": sample, "arm": arm, "species": sp, "estimate": est.get(sp, 0.0),
                             "expected": truth.get(sp), "in_truth": sp in truth})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "real_species.tsv")
    return t


# ─── Credible intervals ─────────────────────────────────────────────────────

def _auc(miss, score):
    """P(score of a miss > score of a hit), ties counted half."""
    from scipy.stats import mannwhitneyu
    miss = np.asarray(miss, bool)
    score = np.asarray(score, float)
    if miss.all() or not miss.any():
        return np.nan
    return mannwhitneyu(score[miss], score[~miss]).statistic / (miss.sum() * (~miss).sum())


def intervals():
    d = pd.read_csv(RESULTS / "intervals.tsv", sep="\t")
    d = d[d["in_truth"]].copy()
    d["depth"] = d["sample"].str.extract(r"_(\d+)k_")[0].astype(int) * 1000
    d["chemistry"] = d["sample"].str.split("_").str[1]
    d["pair_member"] = d["species"].isin({s for pair in PAIRS for s in pair})
    d["error_over_halfwidth"] = (d["estimate"] - d["truth_design"]).abs() / (d["width"] / 2)
    keys = ["method", "depth", "pair_member"]
    t = d.groupby(keys).agg(n=("covers_design", "size"), coverage_design=("covers_design", "mean"),
                            coverage_realised=("covers_realised", "mean"), median_width=("width", "median"),
                            median_error_over_halfwidth=("error_over_halfwidth", "median")).reset_index()
    write_tsv(t, DOCS_DATA / "intervals_coverage.tsv")
    d["abundance_bin"] = pd.cut(d["truth_design"], [0, 0.001, 0.01, 0.05, 1])
    d["srf_bin"] = pd.cut(d["shared_read_fraction"], [-0.01, 0.05, 0.1, 0.2, 0.5, 1.0])
    for name, col in (("abundance", "abundance_bin"), ("shared_read_fraction", "srf_bin")):
        x = (d.groupby(["method", "depth", col], observed=True)["covers_design"].agg(["mean", "size"])
             .rename(columns={"mean": "coverage_design", "size": "n"}).reset_index())
        x[col] = x[col].astype(str)
        write_tsv(x, DOCS_DATA / f"intervals_coverage_by_{name}.tsv")
    rows = []
    for (method, depth), x in d.dropna(subset=["shared_read_fraction"]).groupby(["method", "depth"]):
        rows.append({"method": method, "depth": depth, "n": len(x), "misses": int((~x["covers_design"]).sum()),
                     "auc_shared_read_fraction": _auc(~x["covers_design"], x["shared_read_fraction"])})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "intervals_warning_auc.tsv")
    runtime = pd.read_csv(RESULTS / "intervals_runtime.tsv", sep="\t")
    write_tsv(runtime, DOCS_DATA / "intervals_runtime.tsv")

    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.8), sharey=True)
    for ax, pair in zip(axes, [False, True]):
        for method, x in t[t["pair_member"] == pair].groupby("method"):
            ax.plot(x["depth"], x["coverage_design"], marker="o", ms=3, label=method)
        ax.axhline(0.95, color="grey", lw=0.6, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("reads per sample")
        ax.set_title("pair members" if pair else "other species")
    axes[0].set_ylabel("coverage of the 95% interval")
    axes[0].set_ylim(0, 1)
    axes[0].legend(frameon=False, fontsize=6)
    save_figure(fig, "fig_intervals", t)
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", nargs="+", choices=["ablation", "contamination", "pairs", "lod", "clade", "real", "intervals"])
    ap.add_argument("--samples", help="comma-separated (default: per analysis)")
    args = ap.parse_args()
    ofat = sorted(pd.read_csv(RESULTS / "grid_ofat.tsv", sep="\t")["sample"].unique())
    cell = sorted(pd.read_csv(RESULTS / "grid_cell.tsv", sep="\t")["sample"].unique())
    for what in args.what:
        samples = args.samples.split(",") if args.samples else (cell if what == "contamination" else ofat)
        if what == "ablation":
            ablation()
        elif what == "contamination":
            contamination(samples)
        elif what == "pairs":
            pairs(samples)
        elif what == "lod":
            lod()
        elif what == "clade":
            clade()
        elif what == "intervals":
            intervals()
        else:
            real()
        print("done", what, flush=True)


if __name__ == "__main__":
    main()
