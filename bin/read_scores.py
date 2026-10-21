#!/usr/bin/env python3
"""Score each read against its candidate species (COMPUTE_LOGP).

Writes a JSON with each read's candidate species and log-likelihoods
("log_p_rgs": {read: [[tax_id, ...], [score, ...]]}) and the numbers of
assigned and unassigned reads, as likelihood.log_prob_rgs_dict returns them.
"""
import argparse
import json

from likelihood import log_prob_rgs_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sam", required=True, help="SAM file of the sample's alignments")
    parser.add_argument("--cigar-json", required=True, help="error profile from cigar_profile.py")
    parser.add_argument("--output", required=True, help="output JSON (<sample>_logp_data.json)")
    args = parser.parse_args()

    with open(args.cigar_json) as f:
        data = json.load(f)

    log_probs     = data["log_probs"]
    zero_locs     = data["zero_locs"]
    longest_align = data["longest_align"]

    log_p_rgs, unassigned_count, assigned_count = log_prob_rgs_dict(
        args.sam, log_probs, longest_align, zero_locs
    )

    with open(args.output, "w") as f:
        json.dump({
            "log_p_rgs": log_p_rgs,
            "assigned_count": assigned_count,
            "unassigned_count": unassigned_count
        }, f)


if __name__ == "__main__":
    main()
