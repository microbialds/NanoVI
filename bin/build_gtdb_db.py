#!/usr/bin/env python3
"""
build_gtdb_db.py — Build a NanoVI-compatible database from GTDB SSU data.

Works with any GTDB release that provides the SSU FASTA file (r214, r220, r226...).
The script reads the taxonomy directly from GTDB's FASTA headers, so no separate
taxonomy file is required.

KEY DESIGN: One taxid per species.
  All genomes belonging to the same species receive the same integer taxid.
  This ensures that reads from, e.g., Pseudomonas aeruginosa are aggregated
  correctly during inference instead of being split across thousands of
  genome-level identifiers.

GTDB name normalisation:
  GTDB uses underscore + uppercase suffixes (e.g. _A, _B, _E) to resolve
  polyphyletic groups. These are stripped so that 'Listeria monocytogenes'
  and 'Listeria monocytogenes_B' collapse into the same species taxid.

USAGE:
  # Download GTDB SSU file from:
  #   https://data.gtdb.ecogenomic.org/releases/latest/
  #   File: ssu_all_rXXX.fna.gz  (bacteria + archaea combined)

  python3 build_gtdb_db.py \\
      --ssu ssu_all_r226.fna.gz \\
      --db-name db_gtdb_r226 \\
      --output-dir /path/to/output/db

OUTPUT:
  <output-dir>/species_taxid.fasta   — FASTA with header: taxid:db_name:count
  <output-dir>/taxonomy.tsv          — one row per species

After building, index with minimap2:
  minimap2 -k 27 -d gtdb_index.mmi <output-dir>/species_taxid.fasta
"""

import os
import re
import sys
import gzip
import hashlib
import argparse
from collections import defaultdict


# ─── Taxonomy ranks expected by NanoVI ───────────────────────────────────────
RANKS_OUT = ["tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom"]

# Mapping from GTDB prefix to rank name used by NanoVI
GTDB_PREFIX_TO_RANK = {
    "d__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}


# ─── Parsing ─────────────────────────────────────────────────────────────────

def normalize_gtdb_name(name: str) -> str:
    """
    Remove GTDB-specific clade suffixes from taxon names.

    GTDB uses underscore + uppercase letters to resolve polyphyletic groups,
    e.g. 'Pseudomonas_E', 'Listeria monocytogenes_B', 'Methanobrevibacter_A'.
    These suffixes have no meaning in standard nomenclature and must be stripped
    so that genomes from the same biological species share the same taxid.

    Examples:
      'Pseudomonas_E aeruginosa'     → 'Pseudomonas aeruginosa'
      'Listeria monocytogenes_B'     → 'Listeria monocytogenes'
      'Methanobrevibacter_A smithii' → 'Methanobrevibacter smithii'
      'Escherichia coli'             → 'Escherichia coli'   (unchanged)
    """
    return re.sub(r'_[A-Z]+', '', name).strip()


def parse_gtdb_lineage(lineage_str: str) -> dict:
    """
    Parse a GTDB lineage string into a rank dictionary.

    Example input:
      'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;
       f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli'

    Returns dict with keys: superkingdom, phylum, class, order, family, genus, species.
    Empty string for ranks not present in the lineage string.
    GTDB-specific suffixes (_A, _B, ...) are stripped from all rank names.
    """
    result = {rank: "" for rank in GTDB_PREFIX_TO_RANK.values()}
    for part in lineage_str.strip().split(";"):
        part = part.strip()
        for prefix, rank in GTDB_PREFIX_TO_RANK.items():
            if part.startswith(prefix):
                result[rank] = normalize_gtdb_name(part[len(prefix):])
                break
    return result


def parse_ssu_header(header: str):
    """
    Parse a GTDB SSU FASTA header line.

    Format:
      >genome_id~contig_id d__Bacteria;p__...;s__Species [location=...] [ssu_len=...] [contig_len=...]

    Returns (genome_id: str, lineage: dict) or (None, None) if unparseable.
    """
    header = header.lstrip(">").strip()
    parts = header.split(" ", 1)
    if len(parts) < 2:
        return None, None

    seq_id = parts[0]
    desc = parts[1]
    genome_id = seq_id.split("~")[0]

    # Lineage string: everything before the first '['
    lineage_str = desc.split("[")[0].strip()
    if not lineage_str.startswith("d__"):
        return None, None

    return genome_id, parse_gtdb_lineage(lineage_str)


# ─── File helpers ─────────────────────────────────────────────────────────────

def open_fasta(path: str):
    """Open a FASTA file, transparently handling .gz compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def fasta_records(path: str):
    """Yield (header_line, sequence_str) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


# ─── Main build function ──────────────────────────────────────────────────────

def build_database(ssu_path: str, db_name: str, output_dir: str, min_length: int = 900):
    """
    Build species_taxid.fasta and taxonomy.tsv from a GTDB SSU FASTA file.

    Parameters
    ----------
    ssu_path   : path to GTDB SSU FASTA (.fna or .fna.gz)
    db_name    : label embedded in FASTA headers (e.g. 'db_gtdb_r226')
    output_dir : directory where output files are written
    min_length : minimum 16S sequence length to include [900 bp]
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── Pass 1: collect unique species and their lineages ─────────────────────
    print("Pass 1: Scanning headers to collect unique species...")
    species_to_lineage = {}
    n_total, n_skipped = 0, 0

    for header, _ in fasta_records(ssu_path):
        n_total += 1
        _, lineage = parse_ssu_header(header)
        if lineage is None or not lineage.get("species"):
            n_skipped += 1
            continue
        sp = lineage["species"]
        if sp not in species_to_lineage:
            species_to_lineage[sp] = lineage

    n_species = len(species_to_lineage)
    print(f"  Total sequences scanned : {n_total:,}")
    print(f"  Skipped (no taxonomy)   : {n_skipped:,}")
    print(f"  Unique species found    : {n_species:,}")

    # Assign stable sorted integer taxids (alphabetical → reproducible across runs)
    species_sorted = sorted(species_to_lineage.keys())
    species_to_taxid = {sp: i + 1 for i, sp in enumerate(species_sorted)}

    # ── Pass 2: build FASTA, deduplicating identical sequences per species ────
    print("\nPass 2: Building FASTA (deduplicating identical sequences per species)...")
    seen_hashes = defaultdict(set)
    out_records = []
    n_short, n_dup = 0, 0
    seq_counter = 0

    for header, seq in fasta_records(ssu_path):
        _, lineage = parse_ssu_header(header)
        if lineage is None or not lineage.get("species"):
            continue
        sp = lineage["species"]
        if sp not in species_to_taxid:
            continue
        if len(seq) < min_length:
            n_short += 1
            continue

        taxid = species_to_taxid[sp]
        seq_hash = hashlib.md5(seq.upper().encode()).hexdigest()

        if seq_hash in seen_hashes[taxid]:
            n_dup += 1
            continue

        seen_hashes[taxid].add(seq_hash)
        seq_counter += 1
        out_records.append((f">{taxid}:{db_name}:{seq_counter}", seq))

    print(f"  Skipped (too short <{min_length} bp) : {n_short:,}")
    print(f"  Skipped (duplicate within species)   : {n_dup:,}")
    print(f"  Unique sequences written             : {len(out_records):,}")

    # ── Write species_taxid.fasta ──────────────────────────────────────────────
    fasta_out = os.path.join(output_dir, "species_taxid.fasta")
    print(f"\nWriting {fasta_out} ...")
    with open(fasta_out, "w") as fh:
        for rec_header, seq in out_records:
            fh.write(f"{rec_header}\n{seq}\n")
    print(f"  Done. ({len(out_records):,} sequences)")

    # ── Write taxonomy.tsv ────────────────────────────────────────────────────
    tax_out = os.path.join(output_dir, "taxonomy.tsv")
    print(f"\nWriting {tax_out} ...")
    with open(tax_out, "w") as fh:
        fh.write("\t".join(RANKS_OUT) + "\n")
        for sp in species_sorted:
            lin = species_to_lineage[sp]
            row = [
                str(species_to_taxid[sp]),
                sp,
                lin.get("genus", ""),
                lin.get("family", ""),
                lin.get("order", ""),
                lin.get("class", ""),
                lin.get("phylum", ""),
                lin.get("superkingdom", ""),
            ]
            fh.write("\t".join(row) + "\n")
    print(f"  Done. ({n_species:,} species)")

    print("\n" + "=" * 60)
    print("Database build complete!")
    print(f"  Output directory : {output_dir}")
    print(f"  Species          : {n_species:,}")
    print(f"  Sequences        : {len(out_records):,}")
    print(f"  Avg seq/species  : {len(out_records)/n_species:.1f}")
    print("=" * 60)
    print("\nNext step — index with minimap2:")
    print(f"  minimap2 -k 27 -d {output_dir}/gtdb_index.mmi {fasta_out}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--ssu", required=True,
        help="Path to GTDB SSU FASTA (.fna or .fna.gz). "
             "Download from https://data.gtdb.ecogenomic.org/releases/latest/",
    )
    parser.add_argument(
        "--db-name", default="db_gtdb",
        help="Short label embedded in FASTA headers (e.g. db_gtdb_r226). [db_gtdb]",
    )
    parser.add_argument(
        "--output-dir", default="./db",
        help="Directory where species_taxid.fasta and taxonomy.tsv are written. [./db]",
    )
    parser.add_argument(
        "--min-length", type=int, default=900,
        help="Minimum 16S sequence length in bp to include. [900]",
    )
    args = parser.parse_args()

    if not os.path.exists(args.ssu):
        sys.exit(f"ERROR: SSU file not found: {args.ssu}")

    build_database(
        ssu_path=args.ssu,
        db_name=args.db_name,
        output_dir=args.output_dir,
        min_length=args.min_length,
    )


if __name__ == "__main__":
    main()
