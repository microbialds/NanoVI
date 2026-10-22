#!/usr/bin/env python3
"""NanoVI on the real mock communities (Phase 7a, items 2c-d, 4, 6, 8).

Samples: ONT Zymo 2025.09 (MSPlus x 6, GutMock x 4; sup and hac basecalls), the
Emu paper's Zymo (SRR10391201) and gut mock (SRR14326315), and the preprint's
Zymo mock (SRR23636353), each subsampled to 10,000 reads (seqkit shuffle, seed
20261022; the first 10,000 of real/sub20k). Arms, on the pipeline's steps:

  reference      exact, -N 3, primary, Emu score, VB, floor 1e-3 (the defaults)
  N50            -N 50                        reps / c995   other databases
  per-base       per-base score                em           EM estimator
  noprune        pruning off, floor 0          unk0.02      unknown component, δ 0.02
  pid<x>         alignments below x% identity skipped
  q10 / q15      fastplong mean quality filter (--min_quality; M8)
  tools_new      fastplong 0.7.0 and minimap2 2.31 (D16)
  emudb          NanoVI on Emu's default database (the 'both ways' comparison)

Writes <bench>/results/real/<sample>.tsv per sample and real_nanovi.tsv.
"""
import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, write_tsv
from evaluate import load_estimate, real_metrics

N_READS = 10000
REAL = BENCH / "real"
EMU_DB = BENCH / "downloads" / "comparator_db" / "emu_default"

# arm -> (filter options, database, -N, likelihood options, estimate options)
ARMS = {
    "reference": ({}, "exact", 3, {}, {}),
    "N50": ({}, "exact", 50, {}, {}),
    "c995": ({}, "c995", 3, {}, {}),
    "reps": ({}, "reps", 3, {}, {}),
    "per-base": ({}, "exact", 3, {"score": "per-base"}, {}),
    "em": ({}, "exact", 3, {}, {"estimator": "em"}),
    "noprune": ({}, "exact", 3, {}, {"min_abundance": 0.0, "no_pruning": True}),
    "unk0.02": ({}, "exact", 3, {"unknown_delta": 0.02}, {"unknown": True}),
    "pid97": ({}, "exact", 3, {"min_identity": 97}, {}),
    "q10": ({"min_quality": 10}, "exact", 3, {}, {}),
    "q15": ({"min_quality": 15}, "exact", 3, {}, {}),
    "tools_new": ({"version": "0.7.0"}, "exact", 3, {}, {}),
    "emudb": ({}, "emu_default", 3, {}, {}),
    "emudb_unk0.02": ({}, "emu_default", 3, {"unknown_delta": 0.02}, {"unknown": True}),
    "emudb_pid97": ({}, "emu_default", 3, {"min_identity": 97}, {}),
}
ALIGN_VERSION = {"tools_new": "2.31"}


def db_dir(db):
    return EMU_DB if db == "emu_default" else BENCH / "db" / db


def subsample(name):
    src = REAL / "sub20k" / f"{name}.fastq.gz"
    out = REAL / "sub10k" / f"{name}.fastq.gz"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(src, "rt") as fin, gzip.open(str(out) + ".tmp", "wt", compresslevel=1) as fout:
            for i, line in enumerate(fin):
                if i >= 4 * N_READS:
                    break
                fout.write(line)
        Path(str(out) + ".tmp").rename(out)
    return out


def arm_paths(name, arm):
    filt_opts, db, N, lik_opts, est_opts = ARMS[arm]
    run = REAL / "runs" / name
    ftag = "_".join(f"{k}{v}" for k, v in sorted(filt_opts.items())) or "default"
    atag = f"{db}_N{N}" + (f"_mm{ALIGN_VERSION[arm]}" if arm in ALIGN_VERSION else "")
    ltag = "_".join(f"{k}{v}" for k, v in sorted(lik_opts.items())) or "default"
    filt = run / f"filtered_{ftag}.fastq.gz"
    aln = run / "aln" / f"{ftag}_{atag}.bam"
    npz = run / "lik" / f"{ftag}_{atag}_{ltag}.npz"
    est = run / "est" / f"{arm}.json"
    return filt, aln, npz, est


def run_sample(name, arms, align_threads):
    reads = subsample(name)
    rows = []
    for arm in arms:
        filt_opts, db, N, lik_opts, est_opts = ARMS[arm]
        filt, aln, npz, est = arm_paths(name, arm)
        harness.filter_reads(reads, filt, **filt_opts)
        version = ALIGN_VERSION.get(arm, "2.24")
        index = harness.index_path(db, version)
        if not index.exists():
            harness.build_index(db_dir(db), index, threads=align_threads, version=version)
        harness.align(filt, index, aln, N, threads=align_threads, version=version)
        harness.likelihood(aln, npz, **lik_opts)
        harness.estimate(npz, est, **est_opts)
        abund, below = load_estimate(est, db_dir(db))
        counts = harness.read_counts(npz)
        rows.append({"sample": name, "arm": arm, "tool": "nanovi", "database": db, "below_detection": below,
                     **counts, **real_metrics(abund, name, N_READS)})
    out = BENCH / "results" / "real" / f"{name}.tsv"
    return write_tsv(pd.DataFrame(rows), out)


def samples():
    return sorted(p.name.replace(".fastq.gz", "") for p in (REAL / "sub20k").glob("*.fastq.gz"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--samples", help="comma-separated (default: every real sample)")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--align-threads", type=int, default=4)
    args = ap.parse_args()
    names = args.samples.split(",") if args.samples else samples()
    arms = args.arms.split(",")
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(run_sample, n, arms, args.align_threads): n for n in names}
        for f in as_completed(futures):
            try:
                f.result()
                print(f"done {futures[f]}", flush=True)
            except Exception as e:
                print(f"FAILED {futures[f]}: {e!r}", flush=True)
    files = sorted((BENCH / "results" / "real").glob("*.tsv"))
    write_tsv(pd.concat(pd.read_csv(p, sep="\t") for p in files), BENCH / "results" / "real_nanovi.tsv")


if __name__ == "__main__":
    main()
