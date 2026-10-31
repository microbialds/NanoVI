#!/usr/bin/env python3
"""Estimate species abundances from the per-read scores (RUN_VI).

Reads the per-read scores that read_likelihood.py writes (an npz with
"offsets", "tax_ids" and "scores").

Writes a JSON with the reported abundances ("abundance": {tax_id: value}), the
number of reads whose species were all pruned ("below_detection_count") and a
per-species read summary ("read_summary": {tax_id: {"expected_reads",
"shared_reads"}}, plus "alpha" with --estimator vb), from which
write_abundance_table.py builds the shared_read_fraction column and, with its
--credible-intervals, the ci_lower and ci_upper columns (D27). With --taxonomy
(and --binomial-names, as given to write_abundance_table.py) a read counts as
shared only when its responsibility is split between rows of the table, not
among tax IDs reported in one row (F46).

With --read-assignments (which needs --taxonomy) it also writes the per-read
assignments, <sample>_read_assignments.tsv.gz (read_assignments.py, F55, D37).
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd

import em_estimator
from read_assignments import DEFAULT_MIN_PROBABILITY, write_read_assignments
from read_scores import load_read_scores
from taxonomy import species_row_keys
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
    parser.add_argument("--taxonomy", help="taxonomy.tsv of the database: shared reads are counted between the "
                                           "table's rows (default: between tax IDs)")
    parser.add_argument("--binomial-names", action="store_true",
                        help="with --taxonomy: rows as write_abundance_table.py --binomial-names makes them")
    parser.add_argument("--read-assignments",
                        help="also write each read's table rows and probabilities to this file "
                             "(<sample>_read_assignments.tsv.gz; needs --taxonomy)")
    parser.add_argument("--read-assignments-min-probability", type=float, default=DEFAULT_MIN_PROBABILITY,
                        help="with --read-assignments: the lowest probability listed "
                             f"(default {DEFAULT_MIN_PROBABILITY})")
    args = parser.parse_args()
    if args.read_assignments and not args.taxonomy:
        parser.error("--read-assignments needs --taxonomy: the file names the table's rows")
    if not 0 <= args.read_assignments_min_probability <= 1:
        parser.error("--read-assignments-min-probability must be between 0 and 1")

    offsets, tax_ids, scores, read_names, unassigned_names = load_read_scores(args.scores,
                                                                              names=bool(args.read_assignments))
    if args.read_assignments and (read_names is None or unassigned_names is None):
        sys.exit(f"ERROR: {args.scores} has no read names; write it with read_likelihood.py --read-names "
                 "to write the read assignments")

    if len(offsets) <= 1:
        # No read was scored (F45): nothing to estimate, and every read is unassigned
        print(f"WARNING: {args.scores} has no scored read; writing an empty abundance", file=sys.stderr)
        with open(args.output, "w") as f:
            json.dump({"abundance": {}, "below_detection_count": 0, "read_summary": {}}, f)
        if args.read_assignments:
            write_read_assignments(args.read_assignments, [], unassigned_names, np.zeros((0, 1), dtype=np.int32),
                                   np.zeros((0, 1)), [], {}, args.read_assignments_min_probability)
        return

    # Species in order of first appearance through a set, as the JSON hand-off did,
    # so the columns (and the floating-point sums) are in the same order as before
    db_ids = list({int(tid) for tid in tax_ids})

    # The table row of each species, so that reads split among the tax IDs of one row are not shared (F46)
    species_rows = None
    if args.taxonomy:
        taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", dtype=str).set_index("tax_id")
        species_rows = dict(zip(db_ids, species_row_keys(db_ids, taxonomy_df, args.binomial_names)))

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
            species_rows=species_rows,
            responsibilities=bool(args.read_assignments),
        )
        freq_full, freq_thresh, _, below_detection_count, read_summary = result[:5]

    reported = freq_thresh if freq_thresh else freq_full
    if args.estimator == "em":
        surviving = list(freq_full)
        idx, logL = sparse_likelihood_from_csr(offsets, tax_ids, scores, np.array(surviving, dtype=np.int64))
        resp = em_estimator.responsibilities(idx, logL, [freq_full[t] for t in surviving], surviving,
                                             list(reported))
        read_summary = em_estimator.read_summary(idx, logL, [freq_full[t] for t in surviving], surviving,
                                                 list(reported), species_rows, resp=resp)

    if args.read_assignments:
        if args.estimator == "em":
            (idx_reported, phi), reported_ids = resp, list(reported)
        else:
            idx_reported, phi, reported_ids = result[5]
        write_read_assignments(args.read_assignments, read_names, unassigned_names, idx_reported, phi,
                               reported_ids, species_rows, args.read_assignments_min_probability)

    output = {
        "abundance": reported,
        "below_detection_count": below_detection_count,
        "read_summary": read_summary,
    }
    with open(args.output, "w") as f:
        json.dump(output, f)


if __name__ == "__main__":
    main()
