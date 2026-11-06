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
import math
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
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf8", newline="")
    return open(path, encoding="utf8", newline="")


def _positive(value, where):
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{where}: copy number {value!r} is not a number")
    if not (number > 0 and math.isfinite(number)):
        raise ValueError(f"{where}: copy number {value!r} must be a positive finite number")
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
            species = (row["species"] or "").strip()
            if not species or not (row["copy_number"] or "").strip():
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
            if row["rank"] == "species" and (row["median"] or "").strip():
                taxid = (row["taxid"] or "").strip()
                if taxid in values:
                    raise ValueError(f"{path}: species taxid {taxid!r} is listed more than once")
                values[taxid] = (_positive(row["median"], f"{path} taxid {taxid}"),
                                 (row.get("childcount") or "").strip())
    return values


def read_taxonomy_rows(path):
    """taxonomy.tsv as a list of dicts, ranks the file lacks given as ''."""
    with open(path, encoding="utf8", newline="") as fh:
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
    with open(path, "w", encoding="utf8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COPY_NUMBER_COLUMNS)
        for r in rows:
            w.writerow([r["tax_id"], r["species"], format_number(r["copy_number"]), r["level"],
                        r["n_genomes"], r["source"]])
    counts = {lvl: sum(r["level"] == lvl for r in rows) for lvl in LEVELS}
    print(f"  Copy numbers ({COPY_NUMBER_FILE}): " + ", ".join(f"{lvl} {n:,}" for lvl, n in counts.items()))
    return counts


def read_copy_number_source(gtdb_metadata=None, table=None, rrndb_pantaxa=None, species_of=None):
    """(values, source label, key) from the one copy-number source given, or None without one.

    The builders read it before building, so a source that is missing or fails to parse
    stops the build before its slow steps.
    """
    given = [s for s in (gtdb_metadata, table, rrndb_pantaxa) if s]
    if len(given) > 1:
        raise ValueError("give one copy-number source: GTDB metadata, a table or rrnDB pan-taxa")
    if not given:
        return None
    if gtdb_metadata:
        return (gtdb_metadata_values(gtdb_metadata, species_of),
                "gtdb_complete_genomes:" + ",".join(os.path.basename(p) for p in gtdb_metadata), "species")
    if table:
        return table_values(table), "table:" + os.path.basename(table), "species"
    return rrndb_pantaxa_values(rrndb_pantaxa), "rrndb_pantaxa:" + os.path.basename(rrndb_pantaxa), "tax_id"


def build_copy_numbers(out_dir, taxonomy_path, gtdb_metadata=None, table=None, rrndb_pantaxa=None,
                       species_of=None, source=None):
    """Write <out_dir>/copy_numbers.tsv from one source, or remove an old one when none is given.

    source: read_copy_number_source's result, read before the build; without it the
    source given here is read.
    """
    path = os.path.join(out_dir, COPY_NUMBER_FILE)
    # A database rebuilt without a source, or whose source fails to parse, must not keep an
    # earlier build's values beside the new taxonomy
    if os.path.exists(path):
        os.remove(path)
    if source is None:
        source = read_copy_number_source(gtdb_metadata, table, rrndb_pantaxa, species_of)
    if source is None:
        return None
    values, label, by = source
    rows = fill_copy_numbers(read_taxonomy_rows(taxonomy_path), values, label, by=by)
    write_copy_numbers(rows, path)
    return rows


def read_copy_numbers(path, tax_ids=None):
    """tax_id -> (species, copy number or None, level) from copy_numbers.tsv.

    tax_ids: read and check only these rows (a sample's tax IDs, not the whole database).
    """
    with open(path, encoding="utf8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"tax_id", "species", "copy_number", "level"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: not a copy_numbers.tsv (no {', '.join(sorted(missing))} column)")
        out = {}
        for line, row in enumerate(reader, 2):
            field = lambda k: (row[k] or "").strip()
            if tax_ids is not None and field("tax_id") not in tax_ids:
                continue
            cn, level = field("copy_number"), field("level")
            if cn and level not in LEVELS[:-1]:
                raise ValueError(f"{path} line {line}: level {level!r} is not one of {', '.join(LEVELS[:-1])}")
            # A tax ID listed twice keeps its first row, as taxonomy.label_lineages keeps its first lineage
            out.setdefault(field("tax_id"), (
                field("species"), _positive(cn, f"{path} line {line}") if cn else None, level))
    return out


def tax_id_copy_numbers(freq, copy_numbers, taxonomy_df):
    """tax_id -> (copy number, level) for the tax IDs of freq; None when it has none.

    A tax ID whose species in copy_numbers.tsv differs from taxonomy.tsv's means the
    two files are of different databases: that is an error. A tax ID missing from
    copy_numbers.tsv has no copy number (with a warning).
    """
    names = {}
    wanted = set(map(str, freq))
    if "species" in taxonomy_df.columns:
        species = taxonomy_df["species"]
        for tid, name in species[species.index.astype(str).isin(wanted)].items():
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


def cell_total(cells):
    """The cells of every row with a copy number, summed in the species table's row order.

    The cell table and the group table divide by this same total, so a group of one row
    gets exactly its row's cell share.
    """
    return sum(cells[i][0] for i in sorted(cells) if cells[i])


def write_cell_table(species_df, cells, path, below_detection_count, unassigned_count):
    """Write <sample>_cell-abundance.tsv from the species table and row_cells' values.

    species_df: freq_to_lineage_df's table, indexed by row; cells: row index -> the
    row_cells tuple or None. Rows keep the species table's order.
    """
    import pandas as pd
    from taxonomy import BELOW_DETECTION_LABEL, TAXONOMY_RANKS, UNASSIGNED_LABEL

    taxa = species_df[~species_df["species"].isin([BELOW_DETECTION_LABEL, UNASSIGNED_LABEL])]
    keep = [i for i in taxa.index if cells.get(i)]
    total = cell_total({i: cells.get(i) for i in taxa.index})
    left_out = taxa.loc[[i for i in taxa.index if not cells.get(i)]]
    # Read share of tax IDs without a copy number inside rows that have one: left out of the
    # row's cell share and copy number, though read_abundance is the whole row's
    partial = max(sum(float(taxa.at[i, "abundance"]) - cells[i][1] for i in keep), 0.0)
    table = pd.DataFrame({
        "species": taxa.loc[keep, "species"],
        "abundance": [cells[i][0] / total for i in keep],
        "read_abundance": taxa.loc[keep, "abundance"].astype(float),
        "copy_number": [cells[i][2] for i in keep],
        "copy_number_level": [cells[i][3] for i in keep],
    })
    for rank in TAXONOMY_RANKS[1:]:
        table[rank] = taxa.loc[keep, rank] if rank in taxa.columns else ""
    no_cn = float(left_out["abundance"].astype(float).sum())
    in_rows = (f"; tax IDs without one in rows listed below: read-level abundance {partial:.6g}, left out of "
               f"their row's cell share and copy number" if partial > 1e-12 else "")
    with open(path, "w", encoding="utf8") as fh:
        fh.write("# Cell-level abundance: read-level abundance / 16S copy number (the database's "
                 "copy_numbers.tsv), renormalised over the rows listed; copy_number_level 'species' is "
                 "the species' own value, 'genus' or 'family' the median it was filled from.\n")
        fh.write(f"# Not included: Below detection ({below_detection_count} reads) and Unassigned "
                 f"({unassigned_count} reads), and {len(left_out)} species without a copy number "
                 f"(read-level abundance {no_cn:.6g}){in_rows}.\n")
        table.to_csv(fh, sep="\t", index=False)
    return table
