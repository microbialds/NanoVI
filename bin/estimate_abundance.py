#!/usr/bin/env python3
"""Estimate species abundances from the per-read scores (RUN_VI).

Reads the per-read scores that read_likelihood.py writes (an npz with
"offsets", "tax_ids" and "scores").

Writes a JSON with the reported abundances ("abundance": {tax_id: value}) and
the number of reads whose species were all pruned ("below_detection_count").
With --credible-intervals it also writes "posterior": {tax_id: {"alpha",
"expected_reads", "shared_reads"}}, from which write_abundance_table.py builds
the experimental ci_lower, ci_upper and shared_read_fraction columns.

--no-pruning and --unknown-component are experimental (Phase 7a benchmark, off
by default). With --unknown-component every read also gets the candidate
UNKNOWN_TAX_ID (0) with its floor score from the npz ("unknown_floor", written by
read_likelihood.py --unknown-delta); its abundance is reported under that ID
(docs/design/unknown_component.md).
"""
import argparse
import json

import numpy as np

import em_estimator
from variational_inference_bayesian import bayesian_vi, sparse_likelihood_from_csr

# Tax ID of the experimental unknown-taxon component (database tax IDs start at 1)
UNKNOWN_TAX_ID = 0


def add_unknown_candidate(offsets, tax_ids, scores, floors):
    """CSR arrays with one more candidate per read: UNKNOWN_TAX_ID, scored at the read's floor."""
    n_reads = len(offsets) - 1
    lengths = np.diff(offsets) + 1
    new_offsets = np.zeros(n_reads + 1, dtype=np.int64)
    np.cumsum(lengths, out=new_offsets[1:])
    last = new_offsets[1:] - 1
    is_old = np.ones(int(new_offsets[-1]), dtype=bool)
    is_old[last] = False
    new_tax_ids = np.empty(int(new_offsets[-1]), dtype=np.int64)
    new_scores = np.empty(int(new_offsets[-1]), dtype=np.float64)
    new_tax_ids[is_old], new_scores[is_old] = tax_ids, scores
    new_tax_ids[last], new_scores[last] = UNKNOWN_TAX_ID, floors
    return new_offsets, new_tax_ids, new_scores


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
    parser.add_argument("--credible-intervals", action="store_true",
                        help="VB only, experimental: write the posterior summary for the interval columns")
    parser.add_argument("--no-pruning", action="store_true",
                        help="experimental: prune only below --min-abundance, not below 1/R or 10/R (both estimators)")
    parser.add_argument("--unknown-component", action="store_true",
                        help="experimental: add the unknown-taxon component (needs read_likelihood.py --unknown-delta)")
    args = parser.parse_args()
    if args.credible_intervals and args.estimator != "vb":
        parser.error("--credible-intervals needs --estimator vb")

    with np.load(args.scores) as data:
        offsets, tax_ids, scores = data["offsets"], data["tax_ids"], data["scores"]
        if args.unknown_component:
            if "unknown_floor" not in data:
                parser.error("--unknown-component needs scores written with read_likelihood.py --unknown-delta")
            offsets, tax_ids, scores = add_unknown_candidate(offsets, tax_ids, scores, data["unknown_floor"])

    # Species in order of first appearance through a set, as the JSON hand-off did,
    # so the columns (and the floating-point sums) are in the same order as before
    db_ids = list({int(tid) for tid in tax_ids})

    posterior = None
    if args.estimator == "em":
        tax_id_lists, score_lists = tax_ids.tolist(), scores.tolist()
        log_p_rgs = {
            r: (tax_id_lists[start:end], score_lists[start:end])
            for r, (start, end) in enumerate(zip(offsets[:-1].tolist(), offsets[1:].tolist()))
        }
        freq_full, freq_thresh, below_detection_count = em_estimator.variational_inference_iterations(
            log_p_rgs, db_ids, lli_thresh=args.elbo_threshold, input_threshold=args.min_abundance,
            prune=not args.no_pruning,
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
            posterior_summary=args.credible_intervals,
            prune=not args.no_pruning,
        )
        freq_full, freq_thresh, _, below_detection_count = result[:4]
        if args.credible_intervals:
            posterior = result[4]

    output = {
        "abundance": freq_thresh if freq_thresh else freq_full,
        "below_detection_count": below_detection_count
    }
    if posterior is not None:
        output["posterior"] = posterior
    with open(args.output, "w") as f:
        json.dump(output, f)


if __name__ == "__main__":
    main()
