#!/usr/bin/env python3
"""Estimate species abundances from the per-read scores (RUN_VI).

Reads the per-read scores that read_likelihood.py writes (an npz with
"offsets", "tax_ids" and "scores").

Writes a JSON with the reported abundances ("abundance": {tax_id: value}), the
number of reads whose species were all pruned ("below_detection_count") and a
per-species read summary ("read_summary": {tax_id: {"expected_reads",
"shared_reads"}}, plus "alpha" with --estimator vb), from which
write_abundance_table.py builds the shared_read_fraction column and, with its
--credible-intervals, the ci_lower and ci_upper columns (D27).
"""
import argparse
import json
import sys

import numpy as np

import em_estimator
from variational_inference_bayesian import bayesian_vi, sparse_likelihood_from_csr


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--scores", required=True, help="per-read scores (npz) from read_likelihood.py")
    parser.add_argument("--output", required=True, help="output JSON (<sample>_abundance.json)")
    parser.add_argument("--estimator", choices=["vb", "em"], default="vb",
                        help="'vb': Bayesian variational inference (default); 'em': Emu-style EM")
    parser.add_argument("--min-abundance", type=float, default=0.001,
                        help="lowest abundance reported (default 0.001)")
    parser.add_argument("--tolerance", type=float, default=1e-7,
                        help="VB: stop CAVI when no posterior mean abundance changes by more (default 1e-7)")
    parser.add_argument("--max-iterations", type=int, default=1000,
                        help="VB: most CAVI iterations per pruning round (default 1000)")
    parser.add_argument("--alpha-prior", type=float, default=1.0,
                        help="VB: total concentration of the symmetric Dirichlet prior (default 1.0)")
    parser.add_argument("--elbo-threshold", type=float, default=1e-4,
                        help="EM: stop when the log-likelihood rises by less (default 1e-4)")
    args = parser.parse_args()

    with np.load(args.scores) as data:
        offsets, tax_ids, scores = data["offsets"], data["tax_ids"], data["scores"]

    if len(offsets) <= 1:
        # No read was scored (F45): nothing to estimate, and every read is unassigned
        print(f"WARNING: {args.scores} has no scored read; writing an empty abundance", file=sys.stderr)
        with open(args.output, "w") as f:
            json.dump({"abundance": {}, "below_detection_count": 0, "read_summary": {}}, f)
        return

    # Species in order of first appearance through a set, as the JSON hand-off did,
    # so the columns (and the floating-point sums) are in the same order as before
    db_ids = list({int(tid) for tid in tax_ids})

    if args.estimator == "em":
        tax_id_lists, score_lists = tax_ids.tolist(), scores.tolist()
        log_p_rgs = {
            r: (tax_id_lists[start:end], score_lists[start:end])
            for r, (start, end) in enumerate(zip(offsets[:-1].tolist(), offsets[1:].tolist()))
        }
        freq_full, freq_thresh, below_detection_count = em_estimator.variational_inference_iterations(
            log_p_rgs, db_ids, lli_thresh=args.elbo_threshold, input_threshold=args.min_abundance
        )
    else:
        idx, logL = sparse_likelihood_from_csr(offsets, tax_ids, scores, np.array(db_ids, dtype=np.int64))
        result = bayesian_vi(
            idx, logL, db_ids,
            abundance_threshold=args.min_abundance,
            alpha_prior=args.alpha_prior,
            max_outer_iterations=10,
            max_inner_iterations=args.max_iterations,
            inner_tolerance=args.tolerance,
            posterior_summary=True,
        )
        freq_full, freq_thresh, _, below_detection_count, read_summary = result

    reported = freq_thresh if freq_thresh else freq_full
    if args.estimator == "em":
        surviving = list(freq_full)
        idx, logL = sparse_likelihood_from_csr(offsets, tax_ids, scores, np.array(surviving, dtype=np.int64))
        read_summary = em_estimator.read_summary(idx, logL, [freq_full[t] for t in surviving], surviving,
                                                 list(reported))

    output = {
        "abundance": reported,
        "below_detection_count": below_detection_count,
        "read_summary": read_summary,
    }
    with open(args.output, "w") as f:
        json.dump(output, f)


if __name__ == "__main__":
    main()
