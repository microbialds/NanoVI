#!/usr/bin/env python3
"""Score each read against its candidate species (LIKELIHOOD).

Reads the sample's alignments (SAM or BAM) once: estimates the CIGAR error
profile from the primary mapped alignments, then scores every read against
each species it aligns to (likelihood.score_reads). Writes the per-read scores
as an npz file ("offsets", "tax_ids", "scores": read r's candidates are
tax_ids[offsets[r]:offsets[r + 1]]) for estimate_abundance.py, and the numbers
of assigned and unassigned reads as a small JSON for write_abundance_table.py.

The option after --counts is experimental (Phase 7a benchmark, off by default):
it changes the alignments the error profile comes from.
"""
import argparse
import json

import numpy as np

from likelihood import read_scores_to_arrays, score_reads


def main():
    # No abbreviations: the removed --score would otherwise be taken for --scores
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], allow_abbrev=False)
    parser.add_argument("--alignments", required=True, help="SAM or BAM file of the sample's alignments")
    parser.add_argument("--scores", required=True, help="output npz of per-read scores (<sample>_read_scores.npz)")
    parser.add_argument("--counts", required=True, help="output JSON of read counts (<sample>_read_counts.json)")
    parser.add_argument("--error-profile", choices=["primary", "all"], default="primary",
                        help="experimental: alignments the error profile is estimated from (default primary)")
    args = parser.parse_args()

    log_p_rgs, unassigned_count, assigned_count = score_reads(
        args.alignments, error_profile=args.error_profile)

    offsets, tax_ids, scores = read_scores_to_arrays(log_p_rgs)
    with open(args.scores, "wb") as f:
        np.savez_compressed(f, offsets=offsets, tax_ids=tax_ids, scores=scores)
    with open(args.counts, "w") as f:
        json.dump({"assigned_count": assigned_count, "unassigned_count": unassigned_count}, f)


if __name__ == "__main__":
    main()
