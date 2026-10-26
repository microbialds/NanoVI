#!/usr/bin/env python3
"""Phase 7b confirmation: the shipped pipeline against the Phase 7a reference setting.

`nextflow run main.nf` of the Phase 7b code (fastplong 0.7.0, minimap2 2.31, every
default) on the 80 simulated samples at 1k and 10k reads and the 13 real samples
(ONT 2025.09 sup basecalls, SRR10391201, SRR14326315, SRR23636353; 10,000 reads
each), with two databases built from GTDB r226 by bin/build_gtdb_db.py:

  merged    --merge-gtdb-suffixes: the Phase 7a database (rebuilt byte-identically)
  suffixes  the new default, GTDB's names kept as GTDB writes them (F16, D4)

This script reads both runs' per-sample tables and compares them with the Phase 7a
reference estimates (bin/ at 34d0411, fastplong 0.2.2, minimap2 2.24, same
database as 'merged'), and on the real samples also with Phase 7a's arm for the
new tools (tools_new, D30), with the Phase 7a metric code (evaluate.py). Names in the
'suffixes' run are mapped to their merged names (normalize_gtdb_name) before
scoring, so all three are scored against the same truth; on the simulations the
'suffixes' run is also scored on the GTDB species of the simulated genomes
(raw_* columns).

Writes docs/benchmarks/data/confirm_7b_per_sample.tsv and confirm_7b_summary.tsv.
"""
import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from build_gtdb_db import normalize_gtdb_name
from common import BENCH, DOCS_DATA, write_tsv
from evaluate import all_metrics, compare, load_estimate, load_truth, real_metrics

N_READS = 10000
PLACEHOLDERS = {"Below detection", "Unassigned"}
REAL_SAMPLES = [f"zymo_MSPlus_sup_rep{i}" for i in range(1, 7)] + [f"zymo_GutMock_sup_rep{i}" for i in range(1, 5)] + [
    "emu_zymo_SRR10391201", "preprint_zymo_SRR23636353", "emu_gutmock_SRR14326315"]
KEY_METRICS = ["species_l1", "group_l1", "genus_l1", "species_recall", "species_precision", "species_fp_count",
               "species_fp_abundance"]
REAL_KEY_METRICS = ["species_l1_known", "group_l1_known", "genus_l1_known", "species_recall", "species_precision",
                    "species_fp_count", "species_fp_abundance"]


def sim_samples():
    return sorted(p.name for p in (BENCH / "sims").iterdir() if "_1k_" in p.name or "_10k_" in p.name)


def read_table(out_dir, sample):
    """({species: abundance}, assigned reads, unassigned reads) from a pipeline per-sample table."""
    df = pd.read_csv(Path(out_dir) / "sample_outputs" / f"{sample}_rel-abundance.tsv", sep="\t",
                     keep_default_na=False)
    taxa = df[~df["species"].isin(PLACEHOLDERS)]
    est = {sp: float(a) for sp, a in zip(taxa["species"], taxa["abundance"]) if float(a) > 0}
    counts = df["estimated_counts"].astype(int)
    unassigned = int(counts[df["species"] == "Unassigned"].sum())
    return est, int(counts.sum()) - unassigned, unassigned


def merge_names(est):
    out = {}
    for name, v in est.items():
        out[normalize_gtdb_name(name)] = out.get(normalize_gtdb_name(name), 0.0) + v
    return out


def reference_7a(sample, real):
    """The Phase 7a reference estimate and read counts of a sample."""
    if real:
        run = BENCH / "real" / "runs" / sample
        est_path, counts_path = run / "est" / "reference.json", run / "lik" / "default_exact_N3_default.counts.json"
    else:
        run = BENCH / "runs" / sample
        est_path = run / "est" / "exact_N3_primary_emu_vb_f0.001.json"
        counts_path = run / "lik" / "exact_N3_primary_emu.counts.json"
    est, _ = load_estimate(est_path, BENCH / "db" / "exact")
    counts = json.loads(counts_path.read_text())
    return est, counts["assigned_count"], counts["unassigned_count"]


def tools_arm_7a(sample):
    """Phase 7a's estimate and assigned reads with fastplong 0.7.0 and minimap2 2.31 (arm tools_new, D30)."""
    run = BENCH / "real" / "runs" / sample
    est, _ = load_estimate(run / "est" / "tools_new.json", BENCH / "db" / "exact")
    counts = json.loads((run / "lik" / "version0.7.0_exact_N3_mm2.31_default.counts.json").read_text())
    return est, counts["assigned_count"]


def l1(a, b):
    return sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def raw_truth(sample):
    """Simulated truth under the GTDB species of the simulated genomes (suffixes kept)."""
    copies = pd.read_csv(BENCH / "sims" / "references" / "copies.tsv", sep="\t")
    genomes = set(copies["genome"])
    raw = {}
    with gzip.open(BENCH / "audit" / "ssu_all" / "sequences.tsv.gz", "rt") as fh:
        for row in pd.read_csv(fh, sep="\t", usecols=["genome", "gtdb_species"], chunksize=200000):
            row = row[row["genome"].isin(genomes)]
            raw.update(zip(row["genome"], row["gtdb_species"]))
    species_raw = {sp: raw[g] for sp, g in copies.drop_duplicates("species")[["species", "genome"]].values}
    truth = load_truth(BENCH / "sims" / sample)
    return {species_raw[s]: v for s, v in zip(truth["species"], truth["realised"])}


def group_of(sample, real):
    if not real:
        shape, chem, depth, _ = sample.split("_")
        return f"{depth} {chem} {shape}"
    if sample.startswith("zymo_"):
        return " ".join(sample.split("_")[1:3])
    return sample.split("_")[-1]


def score(est, sample, real):
    return real_metrics(est, sample, N_READS) if real else all_metrics(est, load_truth(BENCH / "sims" / sample))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--merged", required=True, help="output_dir of the run on the merged (Phase 7a) database")
    ap.add_argument("--suffixes", required=True, help="output_dir of the run on the suffix-kept database")
    args = ap.parse_args()

    raw_truths = {}
    rows = []
    for real, samples in ((False, sim_samples()), (True, REAL_SAMPLES)):
        for sample in samples:
            ref_est, ref_assigned, ref_unassigned = reference_7a(sample, real)
            ref_metrics = score(ref_est, sample, real)
            rows.append({"sample": sample, "kind": "real" if real else "sim", "group": group_of(sample, real),
                         "arm": "7a_reference", "assigned": ref_assigned, "unassigned": ref_unassigned,
                         "species_reported": len(ref_est), **ref_metrics})
            for arm, out_dir in (("7b_merged", args.merged), ("7b_suffixes", args.suffixes)):
                est, assigned, unassigned = read_table(out_dir, sample)
                raw = {}
                if arm == "7b_suffixes":
                    if not real:
                        if not raw_truths:
                            raw_truths.update({s: raw_truth(s) for s in sim_samples()})
                        raw = {f"raw_{k}": v for k, v in compare(est, raw_truths[sample], "species").items()}
                    raw["species_split_names"] = len(est) - len(merge_names(est))
                    est = merge_names(est)
                if arm == "7b_merged" and real:
                    tools_est, tools_assigned = tools_arm_7a(sample)
                    raw.update({"same_species_as_7a_tools": set(est) == set(tools_est),
                                "l1_to_7a_tools": l1(est, tools_est),
                                "assigned_minus_7a_tools": assigned - tools_assigned})
                metrics = score(est, sample, real)
                rows.append({
                    "sample": sample, "kind": "real" if real else "sim", "group": group_of(sample, real),
                    "arm": arm, "assigned": assigned, "unassigned": unassigned, "species_reported": len(est),
                    "same_species_as_7a": set(est) == set(ref_est),
                    "l1_to_7a_estimate": l1(est, ref_est),
                    "assigned_minus_7a": assigned - ref_assigned,
                    **{f"delta_{m}": metrics[m] - ref_metrics[m] for m in (REAL_KEY_METRICS if real else KEY_METRICS)},
                    **metrics, **raw,
                })
    per_sample = pd.DataFrame(rows)
    write_tsv(per_sample, DOCS_DATA / "confirm_7b_per_sample.tsv")

    summaries = []
    for (kind, group, arm), d in per_sample.groupby(["kind", "group", "arm"], sort=True):
        metrics = REAL_KEY_METRICS if kind == "real" else KEY_METRICS
        s = {"kind": kind, "group": group, "arm": arm, "samples": len(d),
             **{m: d[m].mean() for m in metrics}, "assigned": d["assigned"].mean()}
        if arm != "7a_reference":
            s.update({
                "same_species_as_7a": int(d["same_species_as_7a"].sum()),
                "max_l1_to_7a_estimate": d["l1_to_7a_estimate"].max(),
                "max_abs_delta_species_l1": d[f"delta_{metrics[0]}"].abs().max(),
                "max_abs_assigned_minus_7a": d["assigned_minus_7a"].abs().max(),
            })
        if arm == "7b_merged" and kind == "real":
            s.update({"same_species_as_7a_tools": int(d["same_species_as_7a_tools"].sum()),
                      "max_l1_to_7a_tools": d["l1_to_7a_tools"].max(),
                      "max_abs_assigned_minus_7a_tools": d["assigned_minus_7a_tools"].abs().max()})
        if arm == "7b_suffixes":
            s["species_split_names"] = d["species_split_names"].mean()
            if kind == "sim":
                s.update({m: d[m].mean() for m in ("raw_species_l1", "raw_species_recall", "raw_species_fp_count")})
        summaries.append(s)
    summary = pd.DataFrame(summaries)
    write_tsv(summary, DOCS_DATA / "confirm_7b_summary.tsv")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(summary.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
