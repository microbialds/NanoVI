#!/usr/bin/env python3
"""Score each read against its candidate species (LIKELIHOOD).

Reads the sample's alignments (SAM or BAM) once: estimates the CIGAR error
profile from the primary mapped alignments, then scores every read against
each species it aligns to (likelihood.score_reads). Writes the per-read scores
as an npz file ("offsets", "tax_ids", "scores": read r's candidates are
tax_ids[offsets[r]:offsets[r + 1]]) for estimate_abundance.py, and the numbers
of assigned and unassigned reads as a small JSON for write_abundance_table.py.

The options after --counts are experimental (Phase 7a benchmark, off by default):
they change the error profile, the score and which alignments are scored, and
--unknown-delta adds the unknown-taxon floor of each read to the npz
("unknown_floor", in read order; docs/design/unknown_component.md).
"""
import argparse
import json

import numpy as np

from likelihood import read_scores_to_arrays, score_reads


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--alignments", required=True, help="SAM or BAM file of the sample's alignments")
    parser.add_argument("--scores", required=True, help="output npz of per-read scores (<sample>_read_scores.npz)")
    parser.add_argument("--counts", required=True, help="output JSON of read counts (<sample>_read_counts.json)")
    parser.add_argument("--error-profile", choices=["primary", "all"], default="primary",
                        help="experimental: alignments the error profile is estimated from (default primary)")
    parser.add_argument("--score", choices=["emu", "per-base"], default="emu",
                        help="experimental: Emu's score (default) or a per-base score with matches")
    parser.add_argument("--min-identity", type=float,
                        help="experimental: skip alignments below this percent identity (Emu's --min-pid)")
    parser.add_argument("--min-aligned-length", type=int,
                        help="experimental: skip alignments with fewer aligned query bases (Emu's --min-align-len)")
    parser.add_argument("--unknown-delta", type=float,
                        help="experimental: write each read's unknown-taxon floor, at the sample's median "
                             "error rate plus this value")
    args = parser.parse_args()

    result = score_reads(args.alignments, error_profile=args.error_profile, score=args.score,
                         min_identity=args.min_identity, min_aligned_length=args.min_aligned_length,
                         unknown_delta=args.unknown_delta)
    log_p_rgs, unassigned_count, assigned_count = result[:3]

    offsets, tax_ids, scores = read_scores_to_arrays(log_p_rgs)
    arrays = {"offsets": offsets, "tax_ids": tax_ids, "scores": scores}
    if args.unknown_delta is not None:
        arrays["unknown_floor"] = result[3]
    with open(args.scores, "wb") as f:
        np.savez_compressed(f, **arrays)
    with open(args.counts, "w") as f:
        json.dump({"assigned_count": assigned_count, "unassigned_count": unassigned_count}, f)


if __name__ == "__main__":
    main()
