#!/usr/bin/env python3
"""
build_gtdb_db.py — Build a NanoVI-compatible database from GTDB SSU data.

Works with any GTDB release that provides the SSU FASTA file (r214, r220, r226...).
The script reads the taxonomy directly from GTDB's FASTA headers, so no separate
taxonomy file is required.

Names are kept as GTDB writes them, clade suffixes included ('Listeria
monocytogenes_B', 'Pseudomonas_E'): GTDB treats such species as distinct (F16,
D4). --merge-gtdb-suffixes strips the suffixes, merging those species into one
tax ID, and writes merged_species.tsv, which lists what was merged.

Besides species_taxid.fasta and taxonomy.tsv, the output directory gets
identical_16s.tsv (species of one genus that share an identical 16S, D24) and
cross_genus_identical.tsv (sequences found in more than one genus, likely
contaminant contigs; D23), from identical_16s.py.

"""

import os
import sys
import hashlib
import argparse
from collections import defaultdict

import identical_16s


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

    GTDB uses underscore + one or two uppercase letters to resolve polyphyletic groups,
    e.g. 'Pseudomonas_E', 'Listeria monocytogenes_B', 'Methanobrevibacter_A'.
    Stripping them gives the names of standard nomenclature, but merges species
    and genera that GTDB keeps apart, so the builder does it only with
    --merge-gtdb-suffixes (F16).

    Examples:
      'Pseudomonas_E aeruginosa'     → 'Pseudomonas aeruginosa'
      'Listeria monocytogenes_B'     → 'Listeria monocytogenes'
      'Methanobrevibacter_A smithii' → 'Methanobrevibacter smithii'
      'Escherichia coli'             → 'Escherichia coli'   (unchanged)

    The same rule groups genera in identical_16s.py.
    """
    return identical_16s.strip_gtdb_suffixes(name)


def merged_lineage(lineage: dict) -> dict:
    """A rank dictionary with GTDB's clade suffixes stripped from every name (--merge-gtdb-suffixes)."""
    return {rank: normalize_gtdb_name(name) for rank, name in lineage.items()}


def parse_gtdb_lineage(lineage_str: str, merge_suffixes: bool = False) -> dict:
    """
    Parse a GTDB lineage string into a rank dictionary.

    Example input:
      'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;
       f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli'

    Returns dict with keys: superkingdom, phylum, class, order, family, genus, species.
    Empty string for ranks not present in the lineage string. GTDB-specific
    suffixes (_A, _B, ...) are kept, or stripped from all rank names with
    merge_suffixes.
    """
    result = {rank: "" for rank in GTDB_PREFIX_TO_RANK.values()}
    for part in lineage_str.strip().split(";"):
        part = part.strip()
        for prefix, rank in GTDB_PREFIX_TO_RANK.items():
            if part.startswith(prefix):
                result[rank] = part[len(prefix):]
                break
    return merged_lineage(result) if merge_suffixes else result


def parse_ssu_header(header: str, merge_suffixes: bool = False):
    """
    Parse a GTDB SSU FASTA header line.

    Format:
      >genome_id~contig_id d__Bacteria;p__...;s__Species [location=...] [ssu_len=...] [contig_len=...]

    Returns (genome_id: str, lineage: dict) or (None, None) if unparseable;
    merge_suffixes as in parse_gtdb_lineage.
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

    return genome_id, parse_gtdb_lineage(lineage_str, merge_suffixes)


# ─── File helpers ─────────────────────────────────────────────────────────────

# One FASTA reader for the builder and identical_16s.py (benchmarks/ import it from here)
open_fasta = identical_16s.open_fasta
fasta_records = identical_16s.fasta_records


# ─── Main build function ──────────────────────────────────────────────────────

def build_database(ssu_path: str, db_name: str, output_dir: str, min_length: int = 900,
                   merge_suffixes: bool = False):
    """
    Build species_taxid.fasta and taxonomy.tsv from a GTDB SSU FASTA file.

    Parameters
    ----------
    ssu_path   : path to GTDB SSU FASTA (.fna or .fna.gz)
    db_name    : label embedded in FASTA headers (e.g. 'db_gtdb_r226')
    output_dir : directory where output files are written
    min_length : minimum 16S sequence length to include [900 bp]
    merge_suffixes : strip GTDB clade suffixes, merging the species and genera
                 they separate, and write merged_species.tsv (F16)
    """
    os.makedirs(output_dir, exist_ok=True)

    # ── Pass 1: collect unique species and their lineages ─────────────────────
    print("Pass 1: Scanning headers to collect unique species...")
    species_to_lineage = {}
    # With merge_suffixes: per merged name, the GTDB species merged into it and their lineages
    merged_sources = defaultdict(dict)
    n_total, n_skipped = 0, 0

    for header, _ in fasta_records(ssu_path):
        n_total += 1
        _, raw = parse_ssu_header(header)
        if raw is None or not raw.get("species"):
            n_skipped += 1
            continue
        lineage = merged_lineage(raw) if merge_suffixes else raw
        sp = lineage["species"]
        if sp not in species_to_lineage:
            species_to_lineage[sp] = lineage
        if merge_suffixes:
            merged_sources[sp].setdefault(raw["species"], raw)

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

    # With merge_suffixes: the GTDB species with a sequence of at least min_length
    with_sequence = set()

    for header, seq in fasta_records(ssu_path):
        _, raw = parse_ssu_header(header)
        lineage = merged_lineage(raw) if raw is not None and merge_suffixes else raw
        if lineage is None or not lineage.get("species"):
            continue
        sp = lineage["species"]
        if sp not in species_to_taxid:
            continue
        if len(seq) < min_length:
            n_short += 1
            continue
        if merge_suffixes:
            with_sequence.add(raw["species"])

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

    # Species that 16S cannot separate, and sequences shared across genera (D23, D24)
    print()
    identical_16s.write_reports(fasta_out, tax_out, os.path.join(output_dir, "identical_16s.tsv"),
                                os.path.join(output_dir, "cross_genus_identical.tsv"))

    if merge_suffixes:
        n_merged = write_merged_species(merged_sources, species_to_taxid,
                                        os.path.join(output_dir, "merged_species.tsv"), with_sequence)
        print(f"  GTDB species merged by suffix stripping: {n_merged:,} names (merged_species.tsv)")

    print("\n" + "=" * 60)
    print("Database build complete!")
    print(f"  Output directory : {output_dir}")
    print(f"  Species          : {n_species:,}")
    print(f"  Sequences        : {len(out_records):,}")
    print(f"  Avg seq/species  : {len(out_records)/n_species:.1f}")
    print("=" * 60)
    print("\nNext step: pass the directory as --db; the pipeline builds the minimap2")
    print("index (with --type and --kmer_size) and keeps it in --index_cache_dir.")
    print("To build it yourself for --db_index, use the same settings, e.g.:")
    print(f"  minimap2 -x map-ont -k 21 -d {output_dir}/species_taxid.map-ont.k21.mmi {fasta_out}")


def write_merged_species(merged_sources, species_to_taxid, path, with_sequence):
    """Write merged_species.tsv: the names that --merge-gtdb-suffixes built from more than one GTDB species.

    Only GTDB species with a sequence in species_taxid.fasta (with_sequence) count:
    one row per name that merges two or more of them, with its tax ID, those GTDB
    species and their genera (suffixes kept), the ranks at which their lineages
    still differ after stripping, and the GTDB species whose lineage taxonomy.tsv
    keeps (the first seen in the input, with or without a sequence). Returns the
    number of rows.
    """
    rows = []
    for sp in sorted(merged_sources):
        first = next(iter(merged_sources[sp]))
        sources = {name: lin for name, lin in merged_sources[sp].items() if name in with_sequence}
        if len(sources) < 2:
            continue
        conflicts = [rank for rank in RANKS_OUT[2:]
                     if len({normalize_gtdb_name(lin[rank]) for lin in sources.values()}) > 1]
        rows.append([
            str(species_to_taxid[sp]), sp, str(len(sources)), ";".join(sorted(sources)),
            ";".join(sorted({lin["genus"] for lin in sources.values()})), ",".join(conflicts), first,
        ])
    with open(path, "w") as fh:
        fh.write("\t".join(["tax_id", "species", "gtdb_species_merged", "gtdb_species", "gtdb_genera",
                            "ranks_in_conflict", "lineage_from"]) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return len(rows)


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
    parser.add_argument(
        "--merge-gtdb-suffixes", action="store_true",
        help="Strip GTDB clade suffixes ('Listeria monocytogenes_B' -> 'Listeria monocytogenes'), "
             "merging the species and genera they separate, and write merged_species.tsv. "
             "Off by default: names are kept as GTDB writes them.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.ssu):
        sys.exit(f"ERROR: SSU file not found: {args.ssu}")

    build_database(
        ssu_path=args.ssu,
        db_name=args.db_name,
        output_dir=args.output_dir,
        min_length=args.min_length,
        merge_suffixes=args.merge_gtdb_suffixes,
    )


if __name__ == "__main__":
    main()
