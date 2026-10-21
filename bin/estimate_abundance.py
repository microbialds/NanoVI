#!/usr/bin/env python3
"""Estimate species abundances from the per-read scores (RUN_VI).

Writes a JSON with the reported abundances ("abundance": {tax_id: value}) and
the number of reads whose species were all pruned ("below_detection_count").
"""
import argparse
import json

from variational_inference_bayesian import variational_inference_iterations


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--logp", required=True, help="per-read scores from read_scores.py")
    parser.add_argument("--output", required=True, help="output JSON (<sample>_abundance.json)")
    args = parser.parse_args()

    with open(args.logp) as f:
        data = json.load(f)

    log_p_rgs = data["log_p_rgs"]

    db_ids = list({tid for read in log_p_rgs.values() for tid in read[0]})

    freq_full, freq_thresh, below_detection_count = variational_inference_iterations(
        log_p_rgs, db_ids, lli_thresh=0.0001, input_threshold=0.001
    )

    with open(args.output, "w") as f:
        json.dump({
            "abundance": freq_thresh if freq_thresh else freq_full,
            "below_detection_count": below_detection_count
        }, f)


if __name__ == "__main__":
    main()
