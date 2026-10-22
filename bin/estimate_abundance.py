#!/usr/bin/env python3
"""Estimate species abundances from the per-read scores (RUN_VI).

Writes a JSON with the reported abundances ("abundance": {tax_id: value}) and
the number of reads whose species were all pruned ("below_detection_count").
With --credible-intervals it also writes "posterior": {tax_id: {"alpha",
"expected_reads", "shared_reads"}}, from which write_abundance_table.py builds
the experimental ci_lower, ci_upper and shared_read_fraction columns.
"""
import argparse
import json

import em_estimator
from variational_inference_bayesian import bayesian_vi_iterations


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--logp", required=True, help="per-read scores from read_scores.py")
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
    args = parser.parse_args()
    if args.credible_intervals and args.estimator != "vb":
        parser.error("--credible-intervals needs --estimator vb")

    with open(args.logp) as f:
        data = json.load(f)

    log_p_rgs = data["log_p_rgs"]

    db_ids = list({tid for read in log_p_rgs.values() for tid in read[0]})

    posterior = None
    if args.estimator == "em":
        freq_full, freq_thresh, below_detection_count = em_estimator.variational_inference_iterations(
            log_p_rgs, db_ids, lli_thresh=args.elbo_threshold, input_threshold=args.min_abundance
        )
    else:
        result = bayesian_vi_iterations(
            log_p_rgs=log_p_rgs,
            db_species_ids=db_ids,
            abundance_threshold=args.min_abundance,
            alpha_prior=args.alpha_prior,
            max_outer_iterations=10,
            max_inner_iterations=args.max_iterations,
            inner_tolerance=args.tolerance,
            posterior_summary=args.credible_intervals,
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
