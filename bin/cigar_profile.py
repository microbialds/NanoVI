#!/usr/bin/env python3
"""Estimate a sample's CIGAR error profile from its alignments (GET_CIGAR).

Writes a JSON with the log-probabilities of insertions, deletions, soft clips
and mismatches over the primary mapped alignments ("log_probs"), the indices of
the operations with no count ("zero_locs") and each read's longest alignment
("longest_align"), as likelihood.get_cigar_op_log_probabilities returns them.
"""
import argparse
import json

from likelihood import get_cigar_op_log_probabilities


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sam", required=True, help="SAM file of the sample's alignments")
    parser.add_argument("--output", required=True, help="output JSON (<sample>_cigar_info.json)")
    args = parser.parse_args()

    log_probs, zero_locs, longest_align = get_cigar_op_log_probabilities(args.sam)

    with open(args.output, "w") as f:
        json.dump({
            "log_probs": log_probs,
            "zero_locs": zero_locs,
            "longest_align": longest_align
        }, f)


if __name__ == "__main__":
    main()
