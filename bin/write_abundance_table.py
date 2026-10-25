#!/usr/bin/env python3
"""Write a sample's abundance table with its lineages (WRITE_OUTPUT).

Joins the abundances from estimate_abundance.py with taxonomy.tsv and writes
<prefix>.tsv: one row per species name (cut to genus and epithet with
--binomial-names, 'Candidatus' names kept whole), then the 'Below detection' and 'Unassigned'
rows, with an estimated_counts column when --counts is given, the
identical_16s_species column when --identical-16s is given (the species of the
same genus that share an identical 16S with each species, from the database's
identical_16s.tsv; empty when the file is ''), shared_read_fraction (from the
read summary of estimate_abundance.py; D27), and with --credible-intervals the
ci_lower and ci_upper columns (sampling intervals of the Bayesian VI's
mean-field posterior, which ignore how shared reads are split; D7, D27).
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
    parser.add_argument("--credible-intervals", action="store_true",
                        help="add ci_lower and ci_upper, 95%% sampling intervals of the mean-field posterior "
                             "(needs abundances from --estimator vb)")
    parser.add_argument("--binomial-names", action="store_true",
                        help="cut species names to genus and epithet ('Candidatus' names keep three words) "
                             "and sum the species that then share a name")
    parser.add_argument("--identical-16s", metavar="FILE",
                        help="identical_16s.tsv of the database: add the identical_16s_species column "
                             "('' writes the column empty)")
    args = parser.parse_args()

    with open(args.abundance) as f:
        abundance_data = json.load(f)
    read_summary = abundance_data.get("read_summary")
    if args.credible_intervals and not (read_summary is not None
                                        and all("alpha" in v for v in read_summary.values())):
        parser.error("--credible-intervals needs abundances from estimate_abundance.py --estimator vb")

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
        posterior=read_summary,
        intervals=args.credible_intervals,
        identical=read_identical(args.identical_16s) if args.identical_16s is not None else None,
        binomial=args.binomial_names,
    )


if __name__ == "__main__":
    main()
