#!/usr/bin/env python3
"""
Generate phyloseq-compatible consolidated tables from NanoVI per-sample outputs.

Reads multiple *_rel-abundance.tsv files produced by the WRITE_OUTPUT process
and generates four phyloseq-compatible TSV files:
  - otu_table_abundance.tsv  (species × samples, relative abundance)
  - otu_table_counts.tsv     (species × samples, estimated counts)
  - tax_table.tsv            (species × taxonomic ranks)
  - sample_metadata.tsv      (sample-level summary statistics)

Inspired by BugBuster's taxonomy_phyloseq.py but tailored to NanoVI output format.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

from taxonomy import TAXONOMY_RANKS as ALL_RANKS, UNASSIGNED_LABEL

TAXONOMY_RANKS = ALL_RANKS[1:]  # ranks above species


def parse_sample_file(filepath):
    """
    Parse a single NanoVI per-sample *_rel-abundance.tsv file.

    Args:
        filepath: Path to the TSV file.

    Returns:
        Tuple of (sample_name, dataframe).
    """
    sample_name = Path(filepath).stem.replace("_rel-abundance", "")
    df = pd.read_csv(filepath, sep="\t", dtype=str)

    # Ensure numeric columns are properly typed
    if "abundance" in df.columns:
        df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce").fillna(0.0)
    if "estimated_counts" in df.columns:
        df["estimated_counts"] = pd.to_numeric(
            df["estimated_counts"], errors="coerce"
        ).fillna(0)

    return sample_name, df


def build_consolidated_tables(sample_data):
    """
    Build phyloseq-compatible tables from parsed per-sample data.

    Rows that share a species name within a sample are summed (F11). Each
    species' taxonomy is its most common lineage across all rows of all
    samples, ties going to the lineage seen first.

    Args:
        sample_data: List of (sample_name, dataframe) tuples.

    Returns:
        Tuple of (otu_abundance, otu_counts, tax_table, sample_metadata) DataFrames.
        otu_counts is None when no sample has an estimated_counts column.
    """
    abundance_dict = {}
    counts_dict = {}
    lineage_counts = {}
    metadata_rows = []
    any_counts = any("estimated_counts" in df.columns for _, df in sample_data)

    for sample_name, df in sample_data:
        has_counts = "estimated_counts" in df.columns

        for _, row in df.iterrows():
            # WRITE_OUTPUT writes unassigned reads as a row named 'Unassigned' (F01);
            # older tables without a species name map to the same row
            species = str(row.get("species", "")).strip()
            if not species or species == "nan" or species.lower() == UNASSIGNED_LABEL.lower():
                species = UNASSIGNED_LABEL

            # Collect abundance, summing rows that share a species name
            abund_val = float(row.get("abundance", 0.0))
            sample_abund = abundance_dict.setdefault(species, {})
            sample_abund[sample_name] = sample_abund.get(sample_name, 0.0) + abund_val

            # Collect counts
            if has_counts:
                count_val = int(float(row.get("estimated_counts", 0)))
                sample_counts = counts_dict.setdefault(species, {})
                sample_counts[sample_name] = sample_counts.get(sample_name, 0) + count_val

            # Tally lineages; the most common one is used for the taxonomy table
            if species != UNASSIGNED_LABEL:
                lineage = []
                for rank in TAXONOMY_RANKS:
                    val = str(row.get(rank, "")).strip()
                    lineage.append(val if val and val != "nan" else "")
                tally = lineage_counts.setdefault(species, {})
                tally[tuple(lineage)] = tally.get(tuple(lineage), 0) + 1

        # Build per-sample metadata
        taxa_detected = sum(
            1
            for sp, vals in abundance_dict.items()
            if vals.get(sample_name, 0) > 0 and sp != UNASSIGNED_LABEL
        )
        total_abundance = sum(
            vals.get(sample_name, 0)
            for sp, vals in abundance_dict.items()
            if sp != UNASSIGNED_LABEL
        )

        if has_counts:
            assigned_counts = sum(
                vals.get(sample_name, 0)
                for sp, vals in counts_dict.items()
                if sp != UNASSIGNED_LABEL
            )
            unassigned_counts = counts_dict.get(UNASSIGNED_LABEL, {}).get(sample_name, 0)
        else:
            assigned_counts = 0
            unassigned_counts = 0

        metadata_rows.append(
            {
                "sample_id": sample_name,
                "total_abundance": round(total_abundance, 8),
                "n_taxa_detected": taxa_detected,
                "assigned_counts": assigned_counts,
                "unassigned_counts": unassigned_counts,
            }
        )

    # Get ordered list of all samples
    all_samples = [name for name, _ in sample_data]

    # Build OTU abundance table (species × samples)
    all_species = sorted(
        [sp for sp in abundance_dict if sp != UNASSIGNED_LABEL]
    ) + [UNASSIGNED_LABEL]

    otu_abundance = pd.DataFrame(
        {
            sample: [abundance_dict.get(sp, {}).get(sample, 0.0) for sp in all_species]
            for sample in all_samples
        },
        index=all_species,
    )
    otu_abundance.index.name = "species"

    # Build OTU counts table (species × samples), only if some sample has counts
    otu_counts = None
    if any_counts:
        no_counts = [name for name, df in sample_data if "estimated_counts" not in df.columns]
        if no_counts:
            print(
                "Warning: no estimated_counts column in "
                f"{', '.join(no_counts)}; their counts are written as 0",
                file=sys.stderr,
            )
        otu_counts = pd.DataFrame(
            {
                sample: [counts_dict.get(sp, {}).get(sample, 0) for sp in all_species]
                for sample in all_samples
            },
            index=all_species,
        )
        otu_counts.index.name = "species"

    # Build taxonomy table (species × ranks); max() keeps the first of tied lineages
    tax_records = {
        sp: dict(zip(TAXONOMY_RANKS, max(tally, key=tally.get)))
        for sp, tally in lineage_counts.items()
    }
    tax_records[UNASSIGNED_LABEL] = {r: "" for r in TAXONOMY_RANKS}
    tax_table = pd.DataFrame.from_dict(tax_records, orient="index")
    tax_table = tax_table.reindex(all_species)
    tax_table.index.name = "species"
    # Ensure column order
    tax_table = tax_table[TAXONOMY_RANKS]

    # Build sample metadata table
    sample_metadata = pd.DataFrame(metadata_rows).set_index("sample_id")

    return otu_abundance, otu_counts, tax_table, sample_metadata


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate phyloseq-compatible tables from NanoVI per-sample outputs"
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        type=Path,
        help="Per-sample *_rel-abundance.tsv files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory)",
    )
    args = parser.parse_args()

    # Validate inputs
    valid_files = [f for f in args.input_files if f.exists()]
    if not valid_files:
        print("Error: No valid input files found", file=sys.stderr)
        sys.exit(1)

    if len(valid_files) < len(args.input_files):
        missing = len(args.input_files) - len(valid_files)
        print(f"Warning: {missing} input file(s) not found", file=sys.stderr)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Parse all sample files
    print(f"Reading {len(valid_files)} per-sample abundance files...")
    sample_data = [parse_sample_file(f) for f in valid_files]

    # Build consolidated tables
    otu_abundance, otu_counts, tax_table, sample_metadata = build_consolidated_tables(
        sample_data
    )

    # Write outputs
    otu_abund_path = args.output_dir / "otu_table_abundance.tsv"
    otu_abundance.to_csv(otu_abund_path, sep="\t")
    print(f"Saved OTU abundance table: {otu_abund_path}")

    if otu_counts is None:
        print(
            "Warning: no input has an estimated_counts column (run with "
            "--keep_counts true); otu_table_counts.tsv not written",
            file=sys.stderr,
        )
    else:
        otu_counts_path = args.output_dir / "otu_table_counts.tsv"
        otu_counts.to_csv(otu_counts_path, sep="\t")
        print(f"Saved OTU counts table: {otu_counts_path}")

    tax_path = args.output_dir / "tax_table.tsv"
    tax_table.to_csv(tax_path, sep="\t")
    print(f"Saved taxonomy table: {tax_path}")

    meta_path = args.output_dir / "sample_metadata.tsv"
    sample_metadata.to_csv(meta_path, sep="\t")
    print(f"Saved sample metadata: {meta_path}")

    print(
        f"\nPhyloseq table generation complete: "
        f"{len(otu_abundance)} taxa × {len(otu_abundance.columns)} samples"
    )


if __name__ == "__main__":
    main()
