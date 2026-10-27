#!/usr/bin/env python3
# Portions of this file are adapted from Emu (https://github.com/treangenlab/emu),
# Copyright (c) 2021 Kristen Curry, Qi Wang, Michael G. Nute, Alona Tyshaieva,
# Qinglong Wu, Tor Savidge, Alexander Dilthey, Todd J. Treangen. MIT License;
# see THIRD_PARTY_LICENSES for the full text.
# Adapted here: run_build_database (Emu's build_database_main).

"""
NanoVI command-line entry points for the steps that are not part of abundance
estimation: build-database, collapse-taxonomy and combine-outputs. The pipeline
calls them; they can also be run on their own. Abundance estimation runs only
through the Nextflow pipeline (decision D1).
"""
import argparse
import os

import identical_16s
from taxonomy import (
    collapse_rank,
    combine_outputs,
    create_nodes_dict,
    create_names_dict,
    create_species_seq2tax_dict,
    create_direct_seq2tax_dict,
    create_unique_seq_dict,
    create_reduced_fasta,
    build_ncbi_taxonomy,
    build_direct_taxonomy
)

from _version import __version__


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)

    subparsers = parser.add_subparsers(dest="subparser_name", help='sub-commands')

    # 1) BUILD-DATABASE SUBCOMMAND
    build_db_parser = subparsers.add_parser("build-database", help="Build custom GTDB database")
    build_db_parser.add_argument('db_name', type=str)
    build_db_parser.add_argument('--sequences', type=str, required=True)
    build_db_parser.add_argument('--seq2tax', type=str, required=True)
    taxonomy_group = build_db_parser.add_mutually_exclusive_group(required=True)
    taxonomy_group.add_argument('--ncbi-taxonomy', type=str)
    taxonomy_group.add_argument('--taxonomy-list', type=str)

    # 2) COLLAPSE-TAXONOMY SUBCOMMAND
    collapse_parser = subparsers.add_parser("collapse-taxonomy", help="Collapse VI output at specified taxonomic rank")
    collapse_parser.add_argument('input_path', type=str)
    collapse_parser.add_argument('rank', type=str)

    # 3) COMBINE-OUTPUTS SUBCOMMAND
    combine_parser = subparsers.add_parser("combine-outputs", help="Combine VI rel abundance outputs to a single table")
    combine_parser.add_argument('dir_path', type=str)
    combine_parser.add_argument('rank', type=str)
    combine_parser.add_argument('--split-tables', action="store_true")
    combine_parser.add_argument('--counts', action="store_true")

    args = parser.parse_args()

    # -------------------
    #   SUBCOMMAND LOGIC
    # -------------------
    if args.subparser_name == "build-database":
        run_build_database(args)

    elif args.subparser_name == "collapse-taxonomy":
        collapse_rank(args.input_path, args.rank)

    elif args.subparser_name == "combine-outputs":
        combine_outputs(args.dir_path, args.rank, args.split_tables, args.counts)


def run_build_database(args):
    if not args.db_name:
        raise ValueError("Please provide a database name.")
    custom_db_path = os.path.join(os.getcwd(), args.db_name)
    if not os.path.exists(custom_db_path):
        os.makedirs(custom_db_path)
    print(f"custom database generating at path: {custom_db_path} ...")

    if args.ncbi_taxonomy:
        dict_names = create_names_dict(os.path.join(args.ncbi_taxonomy, 'names.dmp'))
        dict_nodes = create_nodes_dict(os.path.join(args.ncbi_taxonomy, 'nodes.dmp'))
        seq2tax = create_species_seq2tax_dict(args.seq2tax, dict_nodes)
    else:
        seq2tax = create_direct_seq2tax_dict(args.seq2tax)

    db_unique_ids = set(seq2tax.values())
    dict_fasta = create_unique_seq_dict(args.sequences, seq2tax)
    from Bio import SeqIO
    fasta_records = create_reduced_fasta(dict_fasta, args.db_name)
    SeqIO.write(fasta_records, os.path.join(custom_db_path, 'species_taxid.fasta'), "fasta")

    output_taxonomy_location = os.path.join(custom_db_path, "taxonomy.tsv")
    if args.ncbi_taxonomy:
        # dict_nodes and dict_names were read from the dump above
        build_ncbi_taxonomy(db_unique_ids, dict_nodes, dict_names, output_taxonomy_location)
    else:
        build_direct_taxonomy(db_unique_ids, args.taxonomy_list, output_taxonomy_location)

    # Species that 16S cannot separate, and sequences shared across genera (D23, D24)
    identical_16s.write_reports(
        os.path.join(custom_db_path, 'species_taxid.fasta'), output_taxonomy_location,
        os.path.join(custom_db_path, 'identical_16s.tsv'), os.path.join(custom_db_path, 'cross_genus_identical.tsv'),
    )
    print("Database creation successful")


if __name__ == "__main__":
    main()