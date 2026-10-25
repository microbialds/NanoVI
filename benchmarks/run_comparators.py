#!/usr/bin/env python3
"""Comparator runs on simulated and real samples (Phase 7a, item 4; decision D8).

  Simulated  the 10k samples of seeds 1-3 (12): Kraken2 + Bracken and wf-16s (minimap2
             and kraken2 modes); Emu (about 2 CPU-s per read) on the four seed-1 samples
             of the grid, on its own database and on NanoVI's, and NanoVI (the grid's
             reference setting) on Emu's database (the 'both ways' comparison; real.py
             runs NanoVI on Emu's database for the real mocks)
  Real       the sup basecalls and the three published mocks (13 samples, 10,000 reads):
             the same tools; Emu on NanoVI's database on rep1 of each ONT community
             and the three published mocks

NanoCLUST (real Zymo mocks only, D20) runs one sample at a time on all the cores,
rep1 of each ONT community and the two published Zymo mocks first, then the other
sup replicates until --nanoclust-hours have passed.

Emu and Kraken2 get the reads NanoVI aligns (length 500-2,000); wf-16s gets the
unfiltered reads and applies its own filter (800-2,000). Writes
<bench>/results/comparators.tsv and comparators_species.tsv.gz.
"""
import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import comparators as cmp
import harness
from common import BENCH, write_tsv
from evaluate import all_metrics, load_estimate, load_truth, per_species, real_metrics
from names import load as load_crosswalk, to_gtdb
from real import EMU_DB, REAL, arm_paths, subsample

SIMS = [f"{shape}_{chem}_10k_s{seed:02d}" for seed in (1, 2, 3) for shape in ("even", "log")
        for chem in ("r941", "r1041")]
SIMS_EMU = [s for s in SIMS if s.endswith("_s01")]
REALS = sorted(p.name.replace(".fastq.gz", "") for p in (REAL / "sub20k").glob("*.fastq.gz") if "_hac_" not in p.name)
REALS_BOTH_WAYS = [s for s in REALS if not s.startswith("zymo_") or s.endswith("_rep1")]
NANOCLUST_SAMPLES = (["zymo_MSPlus_sup_rep1", "zymo_GutMock_sup_rep1", "emu_zymo_SRR10391201",
                      "preprint_zymo_SRR23636353"]
                     + [f"zymo_{c}_sup_rep{i}" for i in range(2, 7) for c in ("MSPlus", "GutMock") if c == "MSPlus" or i < 5])
OUT = BENCH / "comparators"


def inputs(sample):
    """(raw reads, length-filtered reads) as NanoVI sees them."""
    if sample in SIMS:
        raw = BENCH / "sims" / sample / "reads.fastq.gz"
        filt = BENCH / "runs" / sample / "filtered.fastq.gz"
    else:
        raw = subsample(sample)
        filt = arm_paths(sample, "reference")[0]
    harness.filter_reads(raw, filt)
    return raw, filt


def jobs():
    """(tool, database, sample) for every per-sample run."""
    out = []
    for s in SIMS + REALS:
        out.append(("kraken2", "ncbi_16s_18s", s))
    for s in SIMS_EMU + REALS:
        out.append(("emu", "emu_default", s))
    for s in SIMS_EMU + REALS_BOTH_WAYS:
        out.append(("emu", "nanovi_exact", s))
    for s in SIMS_EMU:
        out.append(("nanovi", "emu_default", s))
    return out


def run_nanovi(filt, sample, threads):
    """NanoVI's reference setting (-N 3, VB) on Emu's database."""
    run = OUT / "nanovi" / "emu_default"
    aln, npz, est = (run / f"{sample}{ext}" for ext in (".bam", ".npz", ".json"))
    run.mkdir(parents=True, exist_ok=True)
    harness.align(filt, harness.index_path("emu_default"), aln, 3, threads=threads)
    harness.likelihood(aln, npz)
    return harness.estimate(npz, est)


def run_job(tool, db, sample, threads):
    _, filt = inputs(sample)
    if tool == "nanovi":
        return run_nanovi(filt, sample, threads)
    out_dir = OUT / tool / db
    out_dir.mkdir(parents=True, exist_ok=True)
    # each tool names its output after the reads file, so stage it under the sample's name
    reads = OUT / "reads" / f"{sample}.fastq.gz"
    if not reads.exists():
        reads.parent.mkdir(parents=True, exist_ok=True)
        os.link(filt, reads)
    if tool == "emu":
        return cmp.run_emu(reads, db, out_dir, threads)
    return cmp.run_kraken_bracken(reads, out_dir, threads)


def wf16s_dir(samples, name):
    """One sub-directory per sample with its unfiltered reads (hard links: Docker sees them)."""
    d = OUT / "wf16s_input" / name
    for s in samples:
        raw, _ = inputs(s)
        target = d / s / f"{s}.fastq.gz"
        if not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(raw, target)
    return d


def parse_nanovi(path, ncbi_names):
    est = {}
    for k, v in load_estimate(path, EMU_DB)[0].items():
        est[to_gtdb(k, ncbi_names)] = est.get(to_gtdb(k, ncbi_names), 0.0) + v
    return est


def evaluate(ncbi_names):
    rows, species = [], []

    def add(tool, db, sample, est):
        if sample in SIMS:
            truth = load_truth(BENCH / "sims" / sample)
            metrics = all_metrics(est, truth)
            species.extend({"tool": tool, "database": db, "sample": sample, **r} for r in per_species(est, truth))
        else:
            metrics = real_metrics(est, sample, 10000)
            species.extend({"tool": tool, "database": db, "sample": sample, "species": k, "estimate": v}
                           for k, v in est.items())
        rows.append({"tool": tool, "database": db, "sample": sample, **metrics})

    for tool, db, sample in jobs():
        if tool == "nanovi":
            path = OUT / tool / db / f"{sample}.json"
            parse = parse_nanovi
        elif tool == "emu":
            path = OUT / tool / db / f"{sample}_rel-abundance.tsv"
            parse = cmp.parse_emu
        else:
            path = OUT / tool / db / f"{sample}.bracken.tsv"
            parse = cmp.parse_bracken
        if path.exists():
            # NanoVI's database already has GTDB names: mapping them as NCBI names would rename 339
            names = {} if db == "nanovi_exact" else ncbi_names
            add(tool if tool in ("emu", "nanovi") else "kraken2_bracken", db, sample, parse(path, names))
    for classifier in ("minimap2", "kraken2"):
        for group, samples in (("sims", SIMS), ("real", REALS)):
            path = OUT / "wf16s" / f"{classifier}_{group}" / "abundance_table_species.tsv"
            if path.exists():
                for s in samples:
                    add(f"wf16s_{classifier}", "ncbi_16s_18s", s, cmp.parse_wf16s(path, s, ncbi_names))
    for s in NANOCLUST_SAMPLES:
        path = OUT / "nanoclust" / s / f"{s}.nanoclust_out.txt"
        if path.exists():
            add("nanoclust", "ncbi_16S_ribosomal_RNA", s, cmp.parse_nanoclust(path, ncbi_names))
    write_tsv(pd.DataFrame(rows), BENCH / "results" / "comparators.tsv")
    write_tsv(pd.DataFrame(species), BENCH / "results" / "comparators_species.tsv.gz")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=4, help="per Emu or Kraken2 run")
    ap.add_argument("--wf16s-threads", type=int, default=16)
    ap.add_argument("--nanoclust-hours", type=float, default=8.0)
    ap.add_argument("--skip", default="", help="comma-separated: emu, kraken2, nanovi, wf16s, nanoclust")
    ap.add_argument("--evaluate-only", action="store_true")
    args = ap.parse_args()
    skip = set(filter(None, args.skip.split(",")))
    if not args.evaluate_only:
        todo = [j for j in jobs() if j[0] not in skip]
        # Emu on NanoVI's database first: it is the slowest and the 'both ways' comparison needs it
        todo.sort(key=lambda j: (j[1] != "nanovi_exact", j[0] != "emu"))
        with ProcessPoolExecutor(args.workers) as pool:
            futures = {pool.submit(run_job, *j, args.threads): j for j in todo}
            for f in as_completed(futures):
                try:
                    f.result()
                    print("done", *futures[f], flush=True)
                except Exception as e:
                    print("FAILED", *futures[f], repr(e), flush=True)
        if "wf16s" not in skip:
            for classifier in ("minimap2", "kraken2"):
                for group, samples in (("sims", SIMS), ("real", REALS)):
                    try:
                        cmp.run_wf16s(wf16s_dir(samples, group), OUT / "wf16s" / f"{classifier}_{group}",
                                      classifier, args.wf16s_threads)
                        print("done wf16s", classifier, group, flush=True)
                    except Exception as e:
                        print("FAILED wf16s", classifier, group, repr(e), flush=True)
        if "nanoclust" not in skip:
            start = time.time()
            for s in NANOCLUST_SAMPLES:
                if time.time() - start > args.nanoclust_hours * 3600:
                    print("nanoclust: time budget used, stopping before", s, flush=True)
                    break
                try:
                    t0 = time.time()
                    cmp.run_nanoclust(subsample(s), OUT / "nanoclust", args.wf16s_threads)
                    print(f"done nanoclust {s} {time.time() - t0:.0f} s", flush=True)
                except Exception as e:
                    print("FAILED nanoclust", s, repr(e), flush=True)
    evaluate(load_crosswalk())


if __name__ == "__main__":
    main()
