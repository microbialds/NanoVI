#!/usr/bin/env python3
"""Write a sample's abundance table with its lineages (WRITE_OUTPUT).

Joins the abundances from estimate_abundance.py with taxonomy.tsv and writes
<prefix>.tsv: one row per species, then the 'Below detection' and 'Unassigned'
rows, with an estimated_counts column when --counts is given.
"""
import argparse
import json

import pandas as pd

from taxonomy import freq_to_lineage_df


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--abundance", required=True, help="abundance JSON from estimate_abundance.py")
    parser.add_argument("--logp", required=True, help="per-read scores from read_scores.py (for the read counts)")
    parser.add_argument("--taxonomy", required=True, help="taxonomy.tsv of the database")
    parser.add_argument("--output-prefix", required=True, help="output path without .tsv (<sample>_rel-abundance)")
    parser.add_argument("--counts", action="store_true", help="add the estimated_counts column")
    args = parser.parse_args()

    with open(args.abundance) as f:
        abundance_data = json.load(f)

    freq = abundance_data["abundance"]
    below_detection_count = abundance_data["below_detection_count"]

    with open(args.logp) as f:
        logp_data = json.load(f)

    assigned_count = logp_data["assigned_count"]
    unassigned_count = logp_data["unassigned_count"]

    taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", dtype=str).set_index("tax_id")

    freq_to_lineage_df(
        freq,
        tsv_output_path=args.output_prefix,
        taxonomy_df=taxonomy_df,
        assigned_count=assigned_count,
        unassigned_count=unassigned_count,
        counts=args.counts,
        below_detection_count=below_detection_count
    )


if __name__ == "__main__":
    main()
