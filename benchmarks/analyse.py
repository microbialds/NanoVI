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

Every figure (docs/benchmarks/figures/*.png) has its source data in
docs/benchmarks/data/ under the same name.
"""
import argparse
import hashlib
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
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6), sharey=True)
    for ax, depth in zip(axes, [1000, 10000]):
        for (arm, chem), d in data[data["depth"] == depth].groupby(["arm", "chemistry"]):
            label = {"exact_N3_primary_emu_vb_f0.001": "default", "exact_N3_primary_emu_vb_f0_noprune": "no pruning",
                     "exact_N3_primary_per-base_vb_f0.001": "per-base score"}[arm]
            ax.plot(d["reads_bin"].astype(str), d["detected"], marker="o", ms=3,
                    ls="-" if chem == "r1041" else "--", label=f"{label}, {chem}")
        ax.set_title(f"{depth:,} reads")
        ax.set_xlabel("reads of the species")
    axes[0].set_ylabel("fraction detected")
    axes[1].legend(frameon=False, fontsize=6)
    save_figure(fig, "fig_lod", data)
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("what", nargs="+", choices=["ablation", "contamination", "pairs", "lod"])
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
        else:
            lod()
        print("done", what, flush=True)


if __name__ == "__main__":
    main()
