"""
Taxonomy utilities for NanoVI output generation.

Provides functions for converting species frequency vectors into
lineage-annotated abundance tables, collapsing tables to higher
taxonomic ranks, and combining per-sample results into matrices.
"""

# taxonomy.py

import os
import math
import pandas as pd
import pathlib 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Define rank lists consistently
TAXONOMY_RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
RANKS_PRINTOUT = ["tax_id"] + TAXONOMY_RANKS
# In some places we also need the same ranks in order:
RANKS_ORDER = RANKS_PRINTOUT[:]  # copy

def lineage_dict_from_tid(taxid, nodes_dict, names_dict):
    """
    Given a taxid, walk up the tree (using nodes_dict, names_dict) and produce
    a tuple with (tax_id, species, genus, family, order, class, phylum, superkingdom).
    """
    lineage_list = [taxid] + [""] * (len(RANKS_PRINTOUT) - 1)
    while names_dict[taxid] != "root":
        tup = nodes_dict[taxid]
        if tup[1] in RANKS_PRINTOUT:
            idx = RANKS_PRINTOUT.index(tup[1])
            lineage_list[idx] = names_dict[taxid]
        taxid = tup[0]
    return tuple(lineage_list)

def create_nodes_dict(nodes_path):
    """
    Reads NCBI 'nodes.dmp' file and returns a dict:
    { tax_id: (parent_tax_id, rank) }
    """
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    nodes_df = pd.read_csv(nodes_path, sep='|', header=None, dtype=str)[[0, 1, 2]]
    nodes_df.columns = node_headers
    for col in nodes_df.columns:
        nodes_df[col] = nodes_df[col].str.strip()
    return dict(zip(nodes_df['tax_id'], tuple(zip(nodes_df['parent_tax_id'], nodes_df['rank']))))

def create_names_dict(names_path):
    """
    Reads NCBI 'names.dmp' file and returns a dict { tax_id: name_txt } 
    for only 'scientific name' entries.
    """
    name_headers = ['tax_id', 'name_txt', 'name_class']
    names_df = pd.read_csv(names_path, sep='|', index_col=False, header=None, dtype=str).drop([2, 4], axis=1)
    names_df.columns = name_headers
    for col in names_df.columns:
        names_df[col] = names_df[col].str.strip()
    names_df = names_df[names_df["name_class"] == "scientific name"]
    return dict(zip(names_df['tax_id'], names_df['name_txt']))

def get_species_tid(tid, nodes_dict):
    """
    Given an arbitrary taxid, walk up the lineage until we hit a recognized rank (species, genus, etc.).
    Returns that 'species-level' taxid (or whichever rank from TAXONOMY_RANKS is encountered).
    """
    if str(tid) not in nodes_dict:
        raise ValueError(f"Taxid:{tid} not found in nodes file.")
    # Move up the lineage until the rank is in TAXONOMY_RANKS
    while nodes_dict[str(tid)][1] not in TAXONOMY_RANKS:
        tid = nodes_dict[str(tid)][0]
    return tid

def create_species_seq2tax_dict(seq2tax_path, nodes_dict):
    """
    For each sequence ID -> taxid from seq2tax_path,
    find the 'species-level' tid. Returns dict: { sequence_id: species_tid }
    """
    seq2tax_dict, species_id_dict = {}, {}
    with open(seq2tax_path, encoding="utf8") as file:
        for line in file:
            seqid, tid = line.rstrip().split("\t")
            if tid in species_id_dict:
                species_tid = species_id_dict[tid]
            else:
                species_tid = get_species_tid(tid, nodes_dict)
                species_id_dict[tid] = species_tid
            seq2tax_dict[seqid] = species_tid
    return seq2tax_dict

def create_direct_seq2tax_dict(seq2tax_path):
    """
    If we already have a direct mapping of seqid -> final taxid (no lineage walking needed),
    read them from seq2tax_path. Returns dict { seqid: taxid }.
    """
    seq2_taxid = {}
    with open(seq2tax_path, encoding="utf8") as file:
        for line in file:
            seqid, taxid = line.rstrip().split("\t")
            seq2_taxid[seqid] = taxid
    return seq2_taxid

def create_unique_seq_dict(db_fasta_path, seq2tax_dict):
    """
    For each unique sequence (or reverse complement),
    store the taxids and descriptions in a dictionary to remove redundancy.
    Returns: { Seq: { taxid: [list_of_descriptions] } }
    """
    fasta_dict = {}
    for record in SeqIO.parse(db_fasta_path, "fasta"):
        tid = seq2tax_dict.get(record.id)
        if tid:
            # Check if we already have this seq or its revcomp
            if record.seq in fasta_dict:
                if tid in fasta_dict[record.seq]:
                    fasta_dict[record.seq][tid].append(record.description)
                else:
                    fasta_dict[record.seq][tid] = [record.description]
            elif record.seq.reverse_complement() in fasta_dict:
                revcomp = record.seq.reverse_complement()
                if tid in fasta_dict[revcomp]:
                    fasta_dict[revcomp][tid].append(record.description)
                else:
                    fasta_dict[revcomp][tid] = [record.description]
            else:
                # new entry
                fasta_dict[record.seq] = {tid: [record.description]}
    return fasta_dict

def create_reduced_fasta(fasta_dict, db_name):
    """
    Flatten fasta_dict into a list of SeqRecords, giving each unique (seq, tid) pair
    a new ID: e.g. taxid:db_name:count
    """
    records = []
    count = 1
    for seq, tid_dict in fasta_dict.items():
        for taxid, descriptions in tid_dict.items():
            rec_id = f"{taxid}:{db_name}:{count}"
            rec_desc = "; ".join(descriptions)
            records.append(SeqRecord(seq, id=rec_id, description=rec_desc))
            count += 1
    return records

def build_ncbi_taxonomy(unique_tids, nodes_dict, names_dict, filepath):
    """
    Write a 'taxonomy.tsv' with columns: tax_id, species, genus, ..., superkingdom
    for each tid in unique_tids, using NCBI nodes + names dicts.
    """
    with open(filepath, "w", encoding="utf8") as file:
        # Write header
        header_str = "\t".join(RANKS_PRINTOUT) + "\n"
        file.write(header_str)

        for tid in unique_tids:
            lineage_tuple = lineage_dict_from_tid(tid, nodes_dict, names_dict)
            line_str = "\t".join(lineage_tuple) + "\n"
            file.write(line_str)

def build_direct_taxonomy(tid_set, lineage_path, taxonomy_file):
    """
    If we have a pre-existing lineage file (lineage_path),
    filter only the TIDs in tid_set, and write them out to 'taxonomy_file'.
    """
    with open(taxonomy_file, 'w', encoding="utf8") as tax_output_file:
        with open(lineage_path, encoding="utf8") as in_file:
            # Write the first line's columns, but replace the first col name with 'tax_id'
            first_line = in_file.readline()
            tax_output_file.write(f"{RANKS_PRINTOUT[0]}\t")
            # everything after the first tab
            remainder = first_line.split("\t", 1)[1]
            tax_output_file.write(remainder)

            for line in in_file:
                tax_id = line.split("\t", 1)[0]
                if tax_id in tid_set:
                    tax_output_file.write(line)

def unify_species_abundance(df):
    """
    Merges strain-level or subspecies-level entries so that 'species' is only 2 words.
    E.g., 'Escherichia coli O157' -> 'Escherichia coli'.
    Summarizes total 'abundance' by species, then retains genus, family, order, etc.
    """
    if 'species' not in df.columns:
        raise ValueError("Column 'species' not found in DataFrame. Check input format.")

    df['species'] = df['species'].str.split(' ').str[:2].str.join(' ')

    taxonomic_levels = ["genus", "family", "order", "class", "phylum", "superkingdom"]
    for level in taxonomic_levels:
        if level not in df.columns:
            df[level] = "Unknown"

    group_columns = ['species'] + taxonomic_levels
    df_grouped = df.groupby(group_columns, as_index=False).agg({'abundance': 'sum'})
    return df_grouped

def freq_to_lineage_df(freq, tsv_output_path, taxonomy_df, assigned_count, unassigned_count, counts=False):
    """
    Convert frequency dict (taxid -> freq) to a DataFrame, unify taxonomy, and write to .tsv.
    freq: dict { tax_id: freq_val }
    tsv_output_path: base path (no .tsv extension) for final file
    taxonomy_df: DataFrame with columns matching RANKS_PRINTOUT (index=tax_id).
    assigned_count / unassigned_count: number of reads assigned/unassigned
    counts: if True, add 'estimated counts' column
    """
    # Build a DataFrame from freq plus an 'unassigned' row
    results_df = pd.DataFrame(
        zip(list(freq.keys()) + ['unassigned'],
            list(freq.values()) + [0]),
        columns=["tax_id", "abundance"]
    ).set_index('tax_id')

    # Join with the full taxonomy info
    results_df = results_df.join(taxonomy_df, how='left').reset_index()

    # Unify species-level names if 'species' col is present
    if 'species' in results_df.columns:
        results_df = unify_species_abundance(results_df)

    # Optionally add estimated counts
    if counts:
        # Multiply abundance by assigned_count for all but the last row, then add unassigned_count
        # (the last row is 'unassigned')
        assigned_vals = (results_df["abundance"] * assigned_count)[:-1]
        assigned_vals = assigned_vals.astype(int)  # or float if you prefer
        unassigned_series = pd.Series([unassigned_count], index=[len(results_df) - 1])
        counts_series = pd.concat([assigned_vals, unassigned_series], ignore_index=True)
        results_df["estimated counts"] = counts_series

    # Reorder columns
    column_order = [
        "species", "abundance", "genus", "family", "order", "class", "phylum",
        "superkingdom", "estimated counts"
    ]
    for col in column_order:
        if col not in results_df.columns:
            # add it (some might be missing if not "counts")
            results_df[col] = ""

    results_df = results_df[column_order]

    # Write the .tsv
    out_file = f"{tsv_output_path}.tsv"
    results_df.to_csv(out_file, sep="\t", index=False)
    print(f"Final taxonomy report saved as {out_file}")

    return results_df

def collapse_rank(path, rank):
    """
    Collapses an existing .tsv output to a given taxonomic rank by grouping abundance or counts.
    e.g.: python main.py collapse-taxonomy output_rel-abundance.tsv species
    """
    df_emu = pd.read_csv(path, sep='\t')
    if rank not in TAXONOMY_RANKS:
        raise ValueError(f"Specified rank '{rank}' must be in: {TAXONOMY_RANKS}")

    # Determine which columns to keep
    keep_ranks = TAXONOMY_RANKS[TAXONOMY_RANKS.index(rank):]

    # We handle "estimated counts" if present
    has_counts = "estimated counts" in df_emu.columns

    if has_counts:
        # e.g., keep abundance + estimated counts
        df_cols = ["abundance", "estimated counts"] + keep_ranks
        df_emu_copy = df_emu[df_cols].replace({'-': 0})
        df_emu_copy = df_emu_copy.astype({'abundance': float, 'estimated counts': float})
    else:
        # just keep abundance
        df_cols = ["abundance"] + keep_ranks
        df_emu_copy = df_emu[df_cols].replace({'-': 0})
        df_emu_copy = df_emu_copy.astype({'abundance': float})

    # Group by the kept ranks
    df_emu_copy = df_emu_copy.groupby(keep_ranks, dropna=False).sum()

    # Write new file
    base, ext = os.path.splitext(path)
    output_path = f"{base}-{rank}.tsv"
    df_emu_copy.to_csv(output_path, sep='\t')
    print(f"File generated: {output_path}")

def combine_outputs(dir_path, rank, split_files=False, count_table=False):
    """
    Combines multiple 'rel-abundance' .tsv files at the specified rank into a single table.
    Optionally splits the final output into separate "taxonomy" and "abundance" files,
    or includes 'estimated counts' instead of 'abundance'.
    """
    # We'll keep columns from 'rank' onward
    keep_ranks = RANKS_ORDER[RANKS_ORDER.index(rank):]
    df_combined_full = pd.DataFrame(columns=keep_ranks, dtype=str)

    metric = 'abundance'
    if count_table:
        metric = 'estimated counts'

    for file in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file)
        file_extension = pathlib.Path(file).suffix
        if file_extension == '.tsv' and 'rel-abundance' in file:
            name = pathlib.Path(file).stem  # e.g. "sample_rel-abundance"
            name = name.replace('_rel-abundance', '')
            df_sample = pd.read_csv(file_path, sep='\t', dtype=str)

            # Convert metric to numeric
            df_sample[metric] = pd.to_numeric(df_sample[metric], errors='coerce').fillna(0)

            if rank in df_sample.columns and metric in df_sample.columns:
                keep_ranks_sample = [r for r in keep_ranks if r in df_sample.columns]

                # If last row is 'unassigned' in 'tax_id', you might handle that,
                # but it doesn't look like you're merging 'tax_id' here, so skip?

                df_sample_reduced = df_sample[keep_ranks_sample + [metric]].rename(columns={metric: name})
                df_sample_reduced = df_sample_reduced.groupby(keep_ranks_sample, dropna=False).sum().reset_index()
                # Convert new column to numeric
                df_sample_reduced[name] = pd.to_numeric(df_sample_reduced[name], errors='coerce').fillna(0)
                df_combined_full = pd.merge(df_combined_full, df_sample_reduced, how='outer')

    # Sort the combined table by the rank col
    if rank in df_combined_full.columns:
        df_combined_full = df_combined_full.set_index(rank).sort_index().reset_index()

    filename_suffix = ""
    if count_table:
        filename_suffix = "-counts"

    if split_files:
        # Write a "taxonomy" file with the columns from keep_ranks
        tax_out_path = os.path.join(dir_path, f"vi-combined-taxonomy-{rank}.tsv")
        print(f"Combined taxonomy table generated: {tax_out_path}")
        df_combined_full[keep_ranks].to_csv(tax_out_path, sep='\t', index=False)

        # Then write an "abundance" file with the rank col plus sample columns
        keep_ranks_minus_rank = keep_ranks[:]
        keep_ranks_minus_rank.remove(rank)
        abundance_cols = [rank] + [c for c in df_combined_full.columns if c not in keep_ranks_minus_rank]
        abundance_out_path = os.path.join(dir_path, f"vi-combined-abundance-{rank}{filename_suffix}.tsv")
        df_combined_full[abundance_cols].to_csv(abundance_out_path, sep='\t', index=False)
        print(f"Combined abundance table generated: {abundance_out_path}")
    else:
        # Write one combined file
        out_path = os.path.join(dir_path, f"vi-combined-{rank}{filename_suffix}.tsv")
        df_combined_full.to_csv(out_path, sep='\t', index=False)
        print(f"Combined table generated: {out_path}")

    return df_combined_full
