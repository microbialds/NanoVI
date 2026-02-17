#!/usr/bin/env python3
"""
NanoVI standalone CLI entry point.

Provides a command-line interface for running NanoVI subcommands
outside of the Nextflow pipeline. This is primarily used for
development and testing; production use should go through Nextflow.
"""
import sys
import argparse
import os
import pathlib
import pandas as pd

# Import your modules
from alignment import generate_alignments
from variational_inference import (
    get_cigar_op_log_probabilities,
    log_prob_rgs_dict,
    variational_inference_iterations,
    output_read_assignments
)
from taxonomy import (
    freq_to_lineage_df,
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
from utils import (
    timed_function,
    filter_reads,
    gather_fastq_paths,  # <-- Make sure we import gather_fastq_paths
)

__version__ = "1.0"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)

    subparsers = parser.add_subparsers(dest="subparser_name", help='sub-commands')

    # 1) ABUNDANCE SUBCOMMAND
    abundance_parser = subparsers.add_parser("abundance", help="Generate relative abundance estimates")
    abundance_parser.add_argument(
        'input_file', type=str, nargs='+',
        help='One or more paths to FASTQ(.gz) files or directories containing FASTQs'
    )
    abundance_parser.add_argument('--type', '-x', choices=['map-ont', 'map-pb', 'sr'], default='map-ont')
    abundance_parser.add_argument('--min-abundance','-a', type=float, default=0.0000001)
    # Adjust to your actual environment variable or default:
    abundance_parser.add_argument('--db', type=str, default=os.environ.get("GTDB_DATABASE_DIR"))
    abundance_parser.add_argument('--kmer-size', '-k', type=int, default=27)
    abundance_parser.add_argument('--N', '-N', type=int, default=3)
    abundance_parser.add_argument('--K', '-K', type=int, default=4000000000)
    abundance_parser.add_argument('--output-dir', type=str, default="./results")
    abundance_parser.add_argument('--split-prefix', type=str, default="temp")
    abundance_parser.add_argument('--output-basename', type=str)
    abundance_parser.add_argument('--keep-files', action="store_true")
    abundance_parser.add_argument('--keep-counts', action="store_true")
    abundance_parser.add_argument('--keep-read-assignments', action="store_true")
    abundance_parser.add_argument('--output-unclassified', action="store_true")
    abundance_parser.add_argument('--threads', type=int, default=3)
    abundance_parser.add_argument('--min-length', type=int, default=500)
    abundance_parser.add_argument('--max-length', type=int, default=2000)

    # 2) BUILD-DATABASE SUBCOMMAND
    build_db_parser = subparsers.add_parser("build-database", help="Build custom GTDB database")
    build_db_parser.add_argument('db_name', type=str)
    build_db_parser.add_argument('--sequences', type=str, required=True)
    build_db_parser.add_argument('--seq2tax', type=str, required=True)
    taxonomy_group = build_db_parser.add_mutually_exclusive_group(required=True)
    taxonomy_group.add_argument('--ncbi-taxonomy', type=str)
    taxonomy_group.add_argument('--taxonomy-list', type=str)

    # 3) COLLAPSE-TAXONOMY SUBCOMMAND
    collapse_parser = subparsers.add_parser("collapse-taxonomy", help="Collapse VI output at specified taxonomic rank")
    collapse_parser.add_argument('input_path', type=str)
    collapse_parser.add_argument('rank', type=str)

    # 4) COMBINE-OUTPUTS SUBCOMMAND
    combine_parser = subparsers.add_parser("combine-outputs", help="Combine VI rel abundance outputs to a single table")
    combine_parser.add_argument('dir_path', type=str)
    combine_parser.add_argument('rank', type=str)
    combine_parser.add_argument('--split-tables', action="store_true")
    combine_parser.add_argument('--counts', action="store_true")

    args = parser.parse_args()

    # -------------------
    #   SUBCOMMAND LOGIC
    # -------------------
    if args.subparser_name == "abundance":
        run_abundance(args)

    elif args.subparser_name == "build-database":
        run_build_database(args)

    elif args.subparser_name == "collapse-taxonomy":
        collapse_rank(args.input_path, args.rank)

    elif args.subparser_name == "combine-outputs":
        combine_outputs(args.dir_path, args.rank, args.split_tables, args.counts)


def run_abundance(args):
    # Prepare output folder
    if not args.db:
        raise ValueError("Database not specified (use --db or set GTDB_DATABASE_DIR).")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Gather all FASTQ(.gz) from the user-specified paths (files or directories)
    fastq_files = gather_fastq_paths(args.input_file)
    if not fastq_files:
        raise ValueError("No .fastq/.fastq.gz files found in the given paths.")

    print(f"Found {len(fastq_files)} FASTQ files. Processing individually...")
    
    # Read in the taxonomy for inference
    df_taxonomy_path = os.path.join(args.db, "taxonomy.tsv")
    if not os.path.exists(df_taxonomy_path):
        raise FileNotFoundError(f"Expected taxonomy.tsv in {args.db}, not found.")
    df_taxonomy = pd.read_csv(df_taxonomy_path, sep='\t', index_col='tax_id', dtype=str)
    db_species_tids = df_taxonomy.index

    # Loop over each FASTQ/FASTQ.GZ file
    for fastq_path in fastq_files:
        print(f"\n=== Processing sample: {fastq_path} ===")

        # Derive an output basename for this sample
        # e.g. /path/to/SRR23636353.fastq.gz -> "SRR23636353"
        sample_name = pathlib.Path(fastq_path).stem
        # If the file is something like "reads.fastq.gz", .stem might be "reads.fastq"
        # You could strip off a trailing ".fastq" if you want:
        if sample_name.endswith(".fastq"):
            sample_name = sample_name[:-6]  # remove ".fastq"
        elif sample_name.endswith(".fq"):
            sample_name = sample_name[:-3]  # remove ".fq"

        # If the user gave a single --output-basename, incorporate it
        # e.g. "myBaseline_sampleName"
        # Otherwise, just use the sample name
        if args.output_basename:
            out_file = os.path.join(args.output_dir, f"{args.output_basename}_{sample_name}")
        else:
            out_file = os.path.join(args.output_dir, sample_name)

        # FILTER reads
        filtered_file = timed_function(
            f"fastplong filtering for {sample_name}",
            filter_reads,
            fastq_path,                   # input
            f"{out_file}_filtered.fastq", # uncompressed output path (will get .gz)
            args.min_length,
            args.max_length,
            args.output_dir
        )

        # ALIGNMENT
        SAM_FILE = generate_alignments([filtered_file], out_file, args.db, args)

        # INFERENCE
        log_prob_cigar_op, locs_p_cigar_zero, longest_align_dict = timed_function(
            f"Get CIGAR Op Log Probabilities for {sample_name}",
            get_cigar_op_log_probabilities,
            SAM_FILE,
            args.threads
        )
        log_prob_rgs, counts_unassigned, counts_assigned = timed_function(
            f"Log Probability RGS for {sample_name}",
            log_prob_rgs_dict,
            SAM_FILE,
            log_prob_cigar_op,
            longest_align_dict,
            locs_p_cigar_zero
        )

        freq_full, freq_set_thresh, read_dist = timed_function(
            f"Variational Inference Iterations for {sample_name}",
            variational_inference_iterations,
            log_prob_rgs,
            db_species_tids,
            0.0001,
            args.min_abundance
        )

        # WRITE OUTPUT
        timed_function(
            f"Write Output for {sample_name}",
            freq_to_lineage_df,
            freq_full,
            f"{out_file}_rel-abundance",
            df_taxonomy,
            counts_assigned,
            counts_unassigned,
            args.keep_counts
        )

        # Note: keep-read-assignments is not currently functional
        # (variational_inference_iterations returns None for read_dist)
        # if args.keep_read_assignments and read_dist:
        #     output_read_assignments(read_dist, f"{out_file}_read-assignment-distributions")

        # If a thresholded freq was computed
        if freq_set_thresh:
            timed_function(
                f"Write Output Thresholded for {sample_name}",
                freq_to_lineage_df,
                freq_set_thresh,
                f"{out_file}_rel-abundance-threshold-{args.min_abundance}",
                df_taxonomy,
                counts_assigned,
                counts_unassigned,
                args.keep_counts
            )

        # (optional) Output unclassified
        if args.output_unclassified:
            from Bio import SeqIO
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            # Re-implement your "output_unclassified" logic here, if desired

        # Cleanup
        if not args.keep_files:
            if os.path.exists(SAM_FILE):
                os.remove(SAM_FILE)
            # If you want, also remove the .fastq (since you have a .fastq.gz from filtering, etc.)
            # e.g. raw_file_uncompressed = f"{out_file}_filtered.fastq"
            # if os.path.exists(raw_file_uncompressed):
            #     os.remove(raw_file_uncompressed)


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
        dict_nodes = create_nodes_dict(os.path.join(args.ncbi_taxonomy, 'nodes.dmp'))
        dict_names = create_names_dict(os.path.join(args.ncbi_taxonomy, 'names.dmp'))
        build_ncbi_taxonomy(db_unique_ids, dict_nodes, dict_names, output_taxonomy_location)
    else:
        build_direct_taxonomy(db_unique_ids, args.taxonomy_list, output_taxonomy_location)

    print("Database creation successful")


if __name__ == "__main__":
    main()