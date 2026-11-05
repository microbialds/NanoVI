"""
16S copy numbers per species and the cell-level abundance table (F58, D40).

A species' share of 16S reads is its share of 16S genes, not of cells: a species
with ten copies per genome gives about ten times the reads per cell of a species
with one. The database builders write copy_numbers.tsv, one row per tax ID of
taxonomy.tsv:

  tax_id, species   as in taxonomy.tsv
  copy_number       16S copies per genome; empty when the species has none
  level             'species' (its own value), 'genus' or 'family' (the median
                    of the species values of its genus, or else of its family),
                    or 'none'
  n_genomes         the source's genomes (or records) behind a species value
  source            where the species values come from

The species values come from one of three sources:

  GTDB metadata     bac120_metadata / ar53_metadata of the GTDB release: the
                    median ssu_count (16S genes HMMER finds with Barrnap's
                    HMMs, the genome's records in ssu_all) over the species'
                    'Complete Genome' assemblies. Draft assemblies collapse
                    their rRNA operons (GTDB r226: 72% of drafts count 0 or 1),
                    so they are not used. bin/build_gtdb_db.py --gtdb-metadata.
  a table           a TSV with 'species' and 'copy_number' columns, matched on
                    the species name (build_gtdb_db.py and 'main.py
                    build-database' --copy-numbers).
  rrnDB pan-taxa    rrnDB's pan-taxa statistics for NCBI taxonomy, the median
                    of the species row with the tax ID's NCBI taxid ('main.py
                    build-database --rrndb-pantaxa', with --ncbi-taxonomy).

A tax ID without a value takes the median of the species values of its genus
(the genus with the ranks above it, so a genus name reused in two families is two
genera), or else of its family; one with neither, or with a blank genus and
family, has none. Filled values can be fractional (a genus median of 4.5).

write_abundance_table.py --copy-numbers writes <sample>_cell-abundance.tsv: each
row's read-level abundance divided by its copy number and renormalised over the
rows that have one. A row that sums several tax IDs (--binomial-names) takes the
cell-weighted mean of their copy numbers, so its cell abundance is the sum of
theirs; a group of the ambiguity table (--groups-output) is treated the same way.
The Below detection and Unassigned reads, and rows without a copy number, have
no cell abundance and are left out; the file's first lines say how much.
"""

import csv
import gzip
import os
import statistics
import sys

COPY_NUMBER_FILE = "copy_numbers.tsv"
COPY_NUMBER_COLUMNS = ["tax_id", "species", "copy_number", "level", "n_genomes", "source"]
# Ranks a missing value is filled from, in order (D40), and the ranks above that key them
FILL_RANKS = ["genus", "family"]
LINEAGE_RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
LEVELS = ["species"] + FILL_RANKS + ["none"]
CELL_COLUMNS = ["species", "abundance", "read_abundance", "copy_number", "copy_number_level"]
COMPLETE_GENOME = "Complete Genome"


def _open_text(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")


def _positive(value, where):
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{where}: copy number {value!r} is not a number")
    if not number > 0:
        raise ValueError(f"{where}: copy number {value!r} must be positive")
    return number


def gtdb_metadata_values(paths, species_of):
    """Species -> (median ssu_count over its complete genomes, number of those genomes).

    paths: GTDB bac120_metadata / ar53_metadata files (.tsv or .tsv.gz);
    species_of: maps a genome's gtdb_taxonomy string to the species name the
    database uses (build_gtdb_db.py strips the suffixes with --merge-gtdb-suffixes).
    Genomes whose ssu_count is 0 or missing are skipped.
    """
    counts = {}
    needed = ["gtdb_taxonomy", "ssu_count", "ncbi_assembly_level"]
    for path in paths:
        with _open_text(path) as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
            missing = [c for c in needed if c not in header]
            if missing:
                raise ValueError(f"{path}: not a GTDB metadata file (no {', '.join(missing)} column)")
            i_tax, i_ssu, i_level = (header.index(c) for c in needed)
            width = max(i_tax, i_ssu, i_level) + 1
            for line in fh:
                row = line.rstrip("\r\n").split("\t")
                if len(row) < width or row[i_level] != COMPLETE_GENOME:
                    continue
                try:
                    n = float(row[i_ssu])
                except ValueError:
                    continue
                species = species_of(row[i_tax])
                if n > 0 and species:
                    counts.setdefault(species, []).append(n)
    return {sp: (statistics.median(v), len(v)) for sp, v in counts.items()}


def table_values(path):
    """Species -> (copy number, '') from a TSV with 'species' and 'copy_number' columns."""
    values = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"species", "copy_number"} <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: needs 'species' and 'copy_number' columns")
        for line, row in enumerate(reader, 2):
            species = row["species"].strip()
            if not species or not row["copy_number"].strip():
                continue
            if species in values:
                raise ValueError(f"{path}: species {species!r} is listed more than once")
            values[species] = (_positive(row["copy_number"], f"{path} line {line}"), "")
    return values


def rrndb_pantaxa_values(path):
    """NCBI taxid -> (median 16S copies, child count) from rrnDB's pan-taxa statistics (species rows)."""
    values = {}
    with _open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"taxid", "rank", "median"} <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: not an rrnDB pan-taxa file (needs taxid, rank and median columns)")
        for row in reader:
            if row["rank"] == "species" and row["median"].strip():
                values[row["taxid"].strip()] = (_positive(row["median"], f"{path} taxid {row['taxid']}"),
                                                row.get("childcount", "").strip())
    return values


def read_taxonomy_rows(path):
    """taxonomy.tsv as a list of dicts, ranks the file lacks given as ''."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [{"tax_id": r["tax_id"], **{rank: (r.get(rank) or "").strip() for rank in LINEAGE_RANKS}}
            for r in rows]


def _lineage_key(row, rank):
    """The rank's name with the ranks above it, or None when the rank is blank."""
    if not row[rank]:
        return None
    return tuple(row[r] for r in LINEAGE_RANKS[LINEAGE_RANKS.index(rank):])


def fill_copy_numbers(taxonomy_rows, values, source, by="species"):
    """One copy_numbers.tsv row per tax ID: its own value, else the genus or family median.

    values: species name (by='species') or tax ID (by='tax_id') -> (copy number,
    n_genomes). Medians are over the tax IDs that have their own value, each once.
    """
    rows = []
    for t in taxonomy_rows:
        value = values.get(t[by])
        rows.append({"tax_id": t["tax_id"], "species": t["species"],
                     "copy_number": value[0] if value else None, "level": "species" if value else "none",
                     "n_genomes": value[1] if value else "", "source": source})
    for rank in FILL_RANKS:
        # Each species name once, so tax IDs that share a name do not outweigh the others
        known = {}
        for t, r in zip(taxonomy_rows, rows):
            key = _lineage_key(t, rank)
            if r["level"] == "species" and key is not None:
                known.setdefault(key, {}).setdefault(t["species"], r["copy_number"])
        for t, r in zip(taxonomy_rows, rows):
            key = _lineage_key(t, rank)
            if r["level"] == "none" and key in known:
                r["copy_number"], r["level"] = statistics.median(known[key].values()), rank
    return rows


def format_number(x):
    return "" if x is None else f"{x:g}"


def write_copy_numbers(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COPY_NUMBER_COLUMNS)
        for r in rows:
            w.writerow([r["tax_id"], r["species"], format_number(r["copy_number"]), r["level"],
                        r["n_genomes"], r["source"]])
    counts = {lvl: sum(r["level"] == lvl for r in rows) for lvl in LEVELS}
    print(f"  Copy numbers ({COPY_NUMBER_FILE}): " + ", ".join(f"{lvl} {n:,}" for lvl, n in counts.items()))
    return counts


def build_copy_numbers(out_dir, taxonomy_path, gtdb_metadata=None, table=None, rrndb_pantaxa=None,
                       species_of=None):
    """Write <out_dir>/copy_numbers.tsv from one source, or remove an old one when none is given."""
    path = os.path.join(out_dir, COPY_NUMBER_FILE)
    given = [s for s in (gtdb_metadata, table, rrndb_pantaxa) if s]
    if len(given) > 1:
        raise ValueError("give one copy-number source: GTDB metadata, a table or rrnDB pan-taxa")
    if not given:
        # A database rebuilt without a source must not keep an earlier build's values
        if os.path.exists(path):
            os.remove(path)
        return None
    taxonomy_rows = read_taxonomy_rows(taxonomy_path)
    if gtdb_metadata:
        values = gtdb_metadata_values(gtdb_metadata, species_of)
        source, by = "gtdb_complete_genomes:" + ",".join(os.path.basename(p) for p in gtdb_metadata), "species"
    elif table:
        values, source, by = table_values(table), "table:" + os.path.basename(table), "species"
    else:
        values, source, by = rrndb_pantaxa_values(rrndb_pantaxa), "rrndb_pantaxa:" + os.path.basename(rrndb_pantaxa), "tax_id"
    rows = fill_copy_numbers(taxonomy_rows, values, source, by=by)
    write_copy_numbers(rows, path)
    return rows


def read_copy_numbers(path):
    """tax_id -> (species, copy number or None, level) from copy_numbers.tsv."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"tax_id", "species", "copy_number", "level"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: not a copy_numbers.tsv (no {', '.join(sorted(missing))} column)")
        out = {}
        for row in reader:
            cn = row["copy_number"].strip()
            # A tax ID listed twice keeps its first row, as taxonomy.label_lineages keeps its first lineage
            out.setdefault(row["tax_id"].strip(), (
                row["species"].strip(), _positive(cn, f"{path} tax_id {row['tax_id']}") if cn else None,
                row["level"].strip()))
    return out


def tax_id_copy_numbers(freq, copy_numbers, taxonomy_df):
    """tax_id -> (copy number, level) for the tax IDs of freq; None when it has none.

    A tax ID whose species in copy_numbers.tsv differs from taxonomy.tsv's means the
    two files are of different databases: that is an error. A tax ID missing from
    copy_numbers.tsv has no copy number (with a warning).
    """
    names = {}
    if "species" in taxonomy_df.columns:
        for tid, name in taxonomy_df["species"].items():
            # copy_numbers.tsv writes names stripped; a blank name, or one pandas reads as missing
            # ('NA', 'null'), cannot be compared and is not checked (it is an 'Unknown taxid' row)
            if not (name is None or name != name or not str(name).strip()):
                names.setdefault(str(tid), str(name).strip())
    out = {}
    for tid in map(str, freq):
        if tid not in copy_numbers:
            print(f"Warning: tax_id {tid} is not in copy_numbers.tsv; it has no copy number", file=sys.stderr)
            out[tid] = None
            continue
        species, cn, level = copy_numbers[tid]
        if tid in names and names[tid] != species:
            raise ValueError(f"copy_numbers.tsv does not match taxonomy.tsv: tax_id {tid} is "
                             f"{species!r} there and {names[tid]!r} in taxonomy.tsv")
        out[tid] = (cn, level) if cn is not None else None
    return out


def coarsest(levels):
    return max(levels, key=LEVELS.index)


def combine_cells(parts):
    """(cells, reads, copy number, level) of several rows or tax IDs taken together.

    The copy number is reads / cells, the cell-weighted mean of the parts' copy numbers,
    or their common value exactly when they all have the same one.
    """
    cells = sum(p[0] for p in parts)
    reads = sum(p[1] for p in parts)
    numbers = {p[2] for p in parts}
    cn = numbers.pop() if len(numbers) == 1 else (reads / cells if cells > 0 else None)
    return cells, reads, cn, coarsest([p[3] for p in parts])


def row_cells(freq, row_of_id, cn_of_id):
    """Per table row: (cells, reads with a copy number, copy number, coarsest level), or None.

    cells = Σ r_t / c_t over the row's tax IDs that have a copy number c_t (r_t: read-level
    abundance), so the row's copy number, (reads with one) / cells, is the cell-weighted mean
    of its tax IDs' copy numbers. Tax IDs with no abundance do not count; a row none of whose
    tax IDs has a copy number gets None.
    """
    parts = {}
    for tid, r in freq.items():
        tid = str(tid)
        row = row_of_id[tid]
        parts.setdefault(row, [])
        cn = cn_of_id.get(tid)
        if cn is not None and float(r) > 0:
            parts[row].append((float(r) / cn[0], float(r), cn[0], cn[1]))
    return {row: combine_cells(p) if p else None for row, p in parts.items()}


def write_cell_table(species_df, cells, path, below_detection_count, unassigned_count):
    """Write <sample>_cell-abundance.tsv from the species table and row_cells' values.

    species_df: freq_to_lineage_df's table, indexed by row; cells: row index -> the
    row_cells tuple or None. Rows keep the species table's order.
    """
    import pandas as pd
    from taxonomy import BELOW_DETECTION_LABEL, TAXONOMY_RANKS, UNASSIGNED_LABEL

    taxa = species_df[~species_df["species"].isin([BELOW_DETECTION_LABEL, UNASSIGNED_LABEL])]
    keep = [i for i in taxa.index if cells.get(i)]
    total = sum(cells[i][0] for i in keep)
    left_out = taxa.loc[[i for i in taxa.index if not cells.get(i)]]
    # Read share of tax IDs without a copy number inside rows that have one
    partial = sum(float(taxa.at[i, "abundance"]) - cells[i][1] for i in keep)
    table = pd.DataFrame({
        "species": taxa.loc[keep, "species"],
        "abundance": [cells[i][0] / total for i in keep],
        "read_abundance": taxa.loc[keep, "abundance"].astype(float),
        "copy_number": [cells[i][2] for i in keep],
        "copy_number_level": [cells[i][3] for i in keep],
    })
    for rank in TAXONOMY_RANKS[1:]:
        table[rank] = taxa.loc[keep, rank] if rank in taxa.columns else ""
    no_cn = float(left_out["abundance"].astype(float).sum()) + max(partial, 0.0)
    with open(path, "w") as fh:
        fh.write("# Cell-level abundance: read-level abundance / 16S copy number (the database's "
                 "copy_numbers.tsv), renormalised over the rows listed; copy_number_level 'species' is "
                 "the species' own value, 'genus' or 'family' the median it was filled from.\n")
        fh.write(f"# Not included: Below detection ({below_detection_count} reads) and Unassigned "
                 f"({unassigned_count} reads), and {len(left_out)} species without a copy number "
                 f"(read-level abundance {no_cn:.6g}).\n")
        table.to_csv(fh, sep="\t", index=False)
    return table
