#!/usr/bin/env python3
"""Write a sample's abundance table with its lineages (WRITE_OUTPUT).

Joins the abundances from estimate_abundance.py with taxonomy.tsv and writes
<prefix>.tsv: one row per species, then the 'Below detection' and 'Unassigned'
rows, with an estimated_counts column when --counts is given, the
identical_16s_species column when --identical-16s is given (the species of the
same genus that share an identical 16S with each species, from the database's
identical_16s.tsv; empty when the file is ''), and the experimental ci_lower,
ci_upper and shared_read_fraction columns when the abundance JSON holds a
posterior summary (estimate_abundance.py --credible-intervals).
"""
import argparse
import json

import pandas as pd

from taxonomy import IDENTICAL_COLUMN, freq_to_lineage_df


def read_identical(path):
    """tax_id -> the species it shares an identical 16S with, from identical_16s.tsv ('' gives none)."""
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return dict(zip(df["tax_id"], df[IDENTICAL_COLUMN]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--abundance", required=True, help="abundance JSON from estimate_abundance.py")
    parser.add_argument("--read-counts", required=True,
                        help="assigned and unassigned read counts (JSON) from read_likelihood.py")
    parser.add_argument("--taxonomy", required=True, help="taxonomy.tsv of the database")
    parser.add_argument("--output-prefix", required=True, help="output path without .tsv (<sample>_rel-abundance)")
    parser.add_argument("--counts", action="store_true", help="add the estimated_counts column")
    parser.add_argument("--identical-16s", metavar="FILE",
                        help="identical_16s.tsv of the database: add the identical_16s_species column "
                             "('' writes the column empty)")
    args = parser.parse_args()

    with open(args.abundance) as f:
        abundance_data = json.load(f)

    freq = abundance_data["abundance"]
    below_detection_count = abundance_data["below_detection_count"]

    with open(args.read_counts) as f:
        read_counts = json.load(f)

    assigned_count = read_counts["assigned_count"]
    unassigned_count = read_counts["unassigned_count"]

    taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", dtype=str).set_index("tax_id")

    freq_to_lineage_df(
        freq,
        tsv_output_path=args.output_prefix,
        taxonomy_df=taxonomy_df,
        assigned_count=assigned_count,
        unassigned_count=unassigned_count,
        counts=args.counts,
        below_detection_count=below_detection_count,
        posterior=abundance_data.get("posterior"),
        identical=read_identical(args.identical_16s) if args.identical_16s is not None else None,
    )


if __name__ == "__main__":
    main()
