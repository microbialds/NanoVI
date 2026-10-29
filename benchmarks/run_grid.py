#!/usr/bin/env python3
"""Ablation on identical alignments (Phase 7a, item 3: M1, M3, M4, M7, F18 -N).

Factors: database {exact, c995, reps} x -N {3, 10, 50} x error profile {primary,
all} x read score {emu, per-base} x estimator {vb, em} x output floor
(--min-abundance) {0, 1e-4, 1e-3}, plus a diagnostic with the pruning loop off
(floor 0). The reference setting is the pipeline's default: exact, -N 3, primary,
emu, vb, 1e-3.

  --design cell   every combination, on the 10k-read samples of seed 1 (both
                  chemistries and shapes): the only arms that need new alignments
  --design ofat   the reference and each of profile, score, estimator, floor and
                  pruning changed alone, on the reference alignments of every sample

Alignment against the exact database costs about 0.83 CPU-s per read (c995 0.31,
reps 0.16; -N matters little), so database and -N are varied only in the cell.
Each worker takes one sample through alignment, scores, estimates and metrics,
and writes <bench>/results/grid/<design>/<sample>.tsv (and _species.tsv.gz) as
soon as it is done; --collect merges them into grid_<design>.tsv.
"""
import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, write_tsv
from evaluate import all_metrics, load_estimate, load_truth, per_species

DBS = ["exact", "c995", "reps"]
NS = [3, 10, 50]
PROFILES = ["primary", "all"]
SCORES = ["emu", "per-base"]
ESTIMATORS = ["vb", "em"]
FLOORS = [0.0, 1e-4, 1e-3]
REFERENCE = dict(db="exact", N=3, profile="primary", score="emu", estimator="vb", floor=1e-3, prune=True)


def arms(design):
    if design == "cell":
        out = [dict(db=d, N=n, profile=p, score=s, estimator=e, floor=f, prune=True)
               for d, n, p, s, e, f in itertools.product(DBS, NS, PROFILES, SCORES, ESTIMATORS, FLOORS)]
        out += [dict(db=d, N=n, profile=p, score=s, estimator=e, floor=0.0, prune=False)
                for d, n, p, s, e in itertools.product(DBS, NS, PROFILES, SCORES, ESTIMATORS)]
        return out
    out = [dict(REFERENCE)]
    for key, values in (("profile", PROFILES), ("score", SCORES), ("estimator", ESTIMATORS), ("floor", FLOORS)):
        out += [{**REFERENCE, key: v} for v in values if v != REFERENCE[key]]
    out += [{**REFERENCE, "floor": 0.0, "prune": False}, {**REFERENCE, "estimator": "em", "floor": 0.0, "prune": False}]
    return out


def arm_label(a):
    return (f"{a['db']}_N{a['N']}_{a['profile']}_{a['score']}_{a['estimator']}_f{a['floor']:g}"
            + ("" if a["prune"] else "_noprune"))


def paths(sample, a):
    run = BENCH / "runs" / sample
    aln = run / "aln" / f"{a['db']}_N{a['N']}.bam"
    npz = run / "lik" / f"{a['db']}_N{a['N']}_{a['profile']}_{a['score']}.npz"
    est = run / "est" / f"{arm_label(a)}.json"
    return run / "filtered.fastq.gz", aln, npz, est


def sample_fields(sample):
    shape, chem, depth, seed = sample.split("_")
    return {"sample": sample, "shape": shape, "chemistry": chem, "depth": int(depth[:-1]) * 1000, "seed": int(seed[1:])}


def run_sample(sample, design, align_threads):
    """Every arm of one sample, then its metrics. Returns the per-sample result file."""
    out = BENCH / "results" / "grid" / design / f"{sample}.tsv"
    if out.exists():
        return out
    harness.filter_reads(BENCH / "sims" / sample / "reads.fastq.gz", BENCH / "runs" / sample / "filtered.fastq.gz")
    sample_arms = arms(design)
    for a in sample_arms:
        filt, aln, npz, est = paths(sample, a)
        harness.align(filt, harness.index_path(a["db"]), aln, a["N"], threads=align_threads)
        harness.likelihood(aln, npz, error_profile=a["profile"], score=a["score"])
        harness.estimate(npz, est, estimator=a["estimator"], min_abundance=a["floor"], no_pruning=not a["prune"])
    truth = load_truth(BENCH / "sims" / sample)
    rows, sp_rows = [], []
    for a in sample_arms:
        filt, aln, npz, est_path = paths(sample, a)
        est, below = load_estimate(est_path, BENCH / "db" / a["db"])
        base = {**sample_fields(sample), "arm": arm_label(a), **a}
        # Wall times under the shared machine's load: kept for reference, not performance numbers
        times = {f"wall_{k}": float(Path(str(p) + ".time").read_text()) for k, p in
                 (("align", aln), ("likelihood", npz), ("estimate", est_path))}
        rows.append({**base, **all_metrics(est, truth), "below_detection": below, **harness.read_counts(npz),
                     "bam_bytes": aln.stat().st_size, **times})
        sp_rows += [{**base, **r} for r in per_species(est, truth)]
    write_tsv(pd.DataFrame(sp_rows), out.with_name(f"{sample}_species.tsv.gz"))
    return write_tsv(pd.DataFrame(rows), out)


def collect(design):
    d = BENCH / "results" / "grid" / design
    grid = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*.tsv"))])
    species = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*_species.tsv.gz"))])
    write_tsv(grid, BENCH / "results" / f"grid_{design}.tsv")
    write_tsv(species, BENCH / "results" / f"grid_{design}_species.tsv.gz")
    print(f"{len(grid)} rows from {grid['sample'].nunique()} samples -> {BENCH}/results/grid_{design}.tsv")


def default_samples(design):
    sims = sorted(p.name for p in (BENCH / "sims").iterdir() if (p / "reads.fastq.gz").exists())
    if design == "cell":
        return [s for s in sims if "_10k_" in s and s.endswith("_s01")]
    # cheapest first, so that results arrive early: 1k, 10k, then 50k for seeds 1-3
    order = {"1k": 0, "10k": 1, "50k": 2}
    keep = [s for s in sims if not ("_50k_" in s and int(s[-2:]) > 3)]
    return sorted(keep, key=lambda s: (order[s.split("_")[2]], s[-3:], s))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--design", choices=["cell", "ofat"], required=True)
    ap.add_argument("--samples", help="comma-separated sample names (default: by design)")
    ap.add_argument("--workers", type=int, default=4, help="samples in parallel (default 4)")
    ap.add_argument("--align-threads", type=int, default=4, help="minimap2 threads per worker (default 4)")
    ap.add_argument("--collect", action="store_true", help="only merge the per-sample results")
    args = ap.parse_args()
    if not args.collect:
        for a in arms(args.design):
            harness.check_options(dict(error_profile=a["profile"], score=a["score"]),
                                  dict(estimator=a["estimator"], min_abundance=a["floor"], no_pruning=not a["prune"]))
        samples = args.samples.split(",") if args.samples else default_samples(args.design)
        print(f"{len(samples)} samples, {len(arms(args.design))} arms each", flush=True)
        with ProcessPoolExecutor(args.workers) as pool:
            futures = {pool.submit(run_sample, s, args.design, args.align_threads): s for s in samples}
            for f in as_completed(futures):
                try:
                    f.result()
                    print(f"done {futures[f]}", flush=True)
                except Exception as e:  # keep the other samples going; the error is in the step's .log
                    print(f"FAILED {futures[f]}: {e!r}", flush=True)
    collect(args.design)


if __name__ == "__main__":
    main()
