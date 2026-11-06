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
With --groups-output it also writes <prefix>.tsv of the ambiguity groups that
estimate_abundance.py --groups put into the JSON (ambiguity_groups.py, F57, D39).
With --copy-numbers (the database's copy_numbers.tsv) and --cell-output it writes
<prefix>.tsv of cell-level abundances, each row's abundance divided by its 16S copy
number and renormalised (copy_numbers.py, F58, D40), and adds the groups' cell
abundance and copy number to the group table. The species table never changes.
"""
import argparse
import json
import sys

import pandas as pd

from ambiguity_groups import write_groups_table
from copy_numbers import read_copy_numbers, row_cells, tax_id_copy_numbers, write_cell_table
from taxonomy import IDENTICAL_COLUMN, freq_to_lineage_df, table_row_of_ids


def read_identical(path):
    """tax_id -> the species it shares an identical 16S with, from identical_16s.tsv ('' gives none)."""
    if not path:
        return {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if not {"tax_id", IDENTICAL_COLUMN} <= set(df.columns):
        # An empty or truncated table must not fail every sample: the column stays empty
        print(f"Warning: {path} has no tax_id and {IDENTICAL_COLUMN} columns; "
              f"the {IDENTICAL_COLUMN} column is left empty", file=sys.stderr)
        return {}
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
    parser.add_argument("--groups-output", metavar="PREFIX",
                        help="also write the ambiguity groups to PREFIX.tsv (<sample>_groups; needs abundances "
                             "from estimate_abundance.py --groups)")
    parser.add_argument("--copy-numbers", metavar="FILE",
                        help="copy_numbers.tsv of the database (needs --cell-output)")
    parser.add_argument("--cell-output", metavar="PREFIX",
                        help="also write cell-level abundances to PREFIX.tsv (<sample>_cell-abundance; "
                             "needs --copy-numbers); with --groups-output the group table gets their sums")
    args = parser.parse_args()
    if bool(args.copy_numbers) != bool(args.cell_output):
        parser.error("--copy-numbers and --cell-output go together")

    with open(args.abundance) as f:
        abundance_data = json.load(f)
    read_summary = abundance_data.get("read_summary")
    if args.credible_intervals and not (read_summary is not None
                                        and all("alpha" in v for v in read_summary.values())):
        parser.error("--credible-intervals needs abundances from estimate_abundance.py --estimator vb")
    if args.groups_output and "groups" not in abundance_data:
        parser.error("--groups-output needs abundances from estimate_abundance.py --groups")

    freq = abundance_data["abundance"]
    below_detection_count = abundance_data["below_detection_count"]

    with open(args.read_counts) as f:
        read_counts = json.load(f)

    assigned_count = read_counts["assigned_count"]
    unassigned_count = read_counts["unassigned_count"]

    taxonomy_df = pd.read_csv(args.taxonomy, sep="\t", dtype=str).set_index("tax_id")

    species_df = freq_to_lineage_df(
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
    cells = row_of_id = None
    if args.cell_output:
        cn_of_id = tax_id_copy_numbers(freq, read_copy_numbers(args.copy_numbers, set(map(str, freq))), taxonomy_df)
        # The species table row of every tax ID, as the group table finds it
        row_of_id = table_row_of_ids(freq, species_df, taxonomy_df, args.binomial_names)
        cells = row_cells(freq, row_of_id, cn_of_id)
        out_file = f"{args.cell_output}.tsv"
        write_cell_table(species_df, cells, out_file, below_detection_count, unassigned_count)
        print(f"Cell-level abundances saved as {out_file}")
    if args.groups_output:
        out_file = f"{args.groups_output}.tsv"
        write_groups_table(abundance_data["groups"], read_summary, species_df, taxonomy_df, out_file,
                           binomial=args.binomial_names, intervals=args.credible_intervals, cells=cells,
                           row_of_id=row_of_id)
        print(f"Ambiguity groups saved as {out_file}")


if __name__ == "__main__":
    main()
