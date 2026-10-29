#!/usr/bin/env python3
"""Mean-quality filter on the simulated samples (Phase 7a, item 6; M8).

fastplong --mean_qual (the pipeline's --min_quality) against the length-only
reference, on the 10k samples (all seeds). fastplong averages the Phred values
themselves, and Badread's are high for the reads' identity: the R9.4.1-like reads
(about 94% identity) have a median mean of Q23 (5th percentile Q19.5) and the
R10.4.1-like ones Q33 (5th percentile Q31.5), so Q10 and Q12 remove nothing. The
thresholds are therefore set from those distributions, to remove roughly 5-25% and
25-50% of the reads: Q20 and Q22 (R9.4.1), Q32 and Q33 (R10.4.1). The alignments are
the reference alignments of the reads that pass (harness.subset_bam).
Writes <bench>/results/quality_sims.tsv.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

import harness
from common import BENCH, write_tsv
from evaluate import all_metrics, load_estimate, load_truth
from run_grid import sample_fields

QUALITIES = {"r941": [None, 20, 22], "r1041": [None, 32, 33]}


def run_sample(sample):
    run = BENCH / "runs" / sample
    reads = BENCH / "sims" / sample / "reads.fastq.gz"
    ref_aln = run / "aln" / "exact_N3.bam"
    rows = []
    for q in QUALITIES[sample.split("_")[1]]:
        tag = "default" if q is None else f"q{q}"
        filt = run / "filtered.fastq.gz" if q is None else run / f"filtered_{tag}.fastq.gz"
        harness.filter_reads(reads, filt, min_quality=q)
        aln = ref_aln if q is None else run / "aln" / f"exact_N3_{tag}.bam"
        if q is not None:
            harness.subset_bam(ref_aln, filt, aln)
        npz = run / "lik" / ("exact_N3_primary_emu.npz" if q is None else f"exact_N3_{tag}.npz")
        est = run / "est" / ("exact_N3_primary_emu_vb_f0.001.json" if q is None else f"quality_{tag}.json")
        harness.likelihood(aln, npz)
        harness.estimate(npz, est)
        abund, below = load_estimate(est, BENCH / "db" / "exact")
        rows.append({**sample_fields(sample), "filter": tag, "below_detection": below, **harness.read_counts(npz),
                     **all_metrics(abund, load_truth(BENCH / "sims" / sample))})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    samples = sorted(p.name for p in (BENCH / "sims").glob("*_10k_s*"))
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(run_sample, s): s for s in samples}
        for f in as_completed(futures):
            try:
                rows += f.result()
            except Exception as e:
                print("FAILED", futures[f], repr(e), flush=True)
    write_tsv(pd.DataFrame(rows), BENCH / "results" / "quality_sims.tsv")


if __name__ == "__main__":
    main()
