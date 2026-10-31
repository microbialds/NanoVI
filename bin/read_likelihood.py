#!/usr/bin/env python3
"""Score each read against its candidate species (LIKELIHOOD).

Reads the sample's alignments (SAM or BAM) once: estimates the CIGAR error
profile from the primary mapped alignments, then scores every read against
each species it aligns to (likelihood.score_reads). Writes the per-read scores
as an npz file ("offsets", "tax_ids", "scores": read r's candidates are
tax_ids[offsets[r]:offsets[r + 1]]; "read_names" and "unassigned_names" for
the per-read assignments, F55) for estimate_abundance.py, and the numbers of
assigned and unassigned reads as a small JSON for write_abundance_table.py.
"""
import argparse
import json
import sys

from likelihood import score_reads
from read_scores import read_scores_to_arrays, save_read_scores


def main():
    # No abbreviations: the removed --score would otherwise be taken for --scores
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], allow_abbrev=False)
    parser.add_argument("--alignments", required=True, help="SAM or BAM file of the sample's alignments")
    parser.add_argument("--scores", required=True, help="output npz of per-read scores (<sample>_read_scores.npz)")
    parser.add_argument("--counts", required=True, help="output JSON of read counts (<sample>_read_counts.json)")
    args = parser.parse_args()

    log_p_rgs, unassigned_count, assigned_count, unassigned_names = score_reads(args.alignments,
                                                                                unassigned_names=True)
    if not assigned_count:
        # The sample still gets a table, with every read unassigned (F45)
        print(f"WARNING: {args.alignments} has no primary alignment to the database; "
              f"its {unassigned_count} reads are unassigned", file=sys.stderr)

    offsets, tax_ids, scores = read_scores_to_arrays(log_p_rgs)
    with open(args.scores, "wb") as f:
        save_read_scores(f, offsets, tax_ids, scores, read_names=list(log_p_rgs), unassigned_names=unassigned_names)
    with open(args.counts, "w") as f:
        json.dump({"assigned_count": assigned_count, "unassigned_count": unassigned_count}, f)


if __name__ == "__main__":
    main()
