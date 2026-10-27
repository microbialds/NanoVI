#!/usr/bin/env python3
"""Species that share an identical 16S sequence in a database (D23, D24).

Reads a database's species_taxid.fasta (headers 'taxid:db_name:n') and
taxonomy.tsv and writes:

  identical_16s.tsv          one row per species that shares at least one
                             identical sequence with another species of its
                             genus: its group (the species of the genus linked
                             through shared sequences) and the species it
                             shares a sequence with. 16S cannot separate them.
  cross_genus_identical.tsv  one row per sequence found in species of more than
                             one genus. In GTDB these are almost always
                             contaminant contigs in an assembly; such a
                             sequence links only the species of each genus
                             that carry it, never species across genera.

Genera are compared with GTDB's clade suffixes stripped ('Pseudomonas_E' is
'Pseudomonas'), so a GTDB database gives the same groups whether or not it was
built with --merge-gtdb-suffixes. Species with no genus in taxonomy.tsv are
left out of both files. build_gtdb_db.py and 'main.py build-database' write
both files into the database directory; the pipeline builds identical_16s.tsv
for a database that lacks it (IDENTICAL_16S) and names the sharing species in
the identical_16s_species column of each sample's table.
"""
import argparse
import csv
import gzip
import hashlib
import re
from collections import defaultdict

# GTDB clade suffix ('Pseudomonas_E', 'Listeria monocytogenes_B'), at the end of a word
# only, so other databases' names ('Lachnospiraceae_NK4A136_group') stay whole;
# build_gtdb_db.py strips it with strip_gtdb_suffixes too
GTDB_SUFFIX = re.compile(r"_[A-Z]+(?=\s|$)")

GROUP_COLUMNS = ["tax_id", "species", "genus", "group", "group_size", "identical_16s_species"]
CROSS_GENUS_COLUMNS = ["md5", "species_count", "sequence_ids", "tax_ids", "species", "genera", "families", "phyla"]


def strip_gtdb_suffixes(name):
    """A taxon name without GTDB's clade suffixes ('Pseudomonas_E aeruginosa' -> 'Pseudomonas aeruginosa')."""
    return GTDB_SUFFIX.sub("", name).strip()


def open_fasta(path):
    """Open a FASTA file, transparently handling .gz compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def fasta_records(path):
    """Yield (header line, sequence) from a FASTA file, gzipped or not (build_gtdb_db.py uses it too)."""
    header, parts = None, []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line, []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def read_taxonomy(path):
    """tax_id -> {species, genus, family, phylum} from taxonomy.tsv."""
    with open(path, newline="", encoding="utf8") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        return {row["tax_id"]: {rank: (row.get(rank) or "").strip() for rank in ("species", "genus", "family", "phylum")}
                for row in rows}


def sequence_carriers(fasta_path):
    """md5 of each distinct sequence -> (tax IDs carrying it, record IDs), in file order."""
    carriers = defaultdict(lambda: ({}, []))
    for header, seq in fasta_records(fasta_path):
        seq_id = header[1:].split()[0]
        tax_ids, seq_ids = carriers[hashlib.md5(seq.upper().encode()).hexdigest()]
        tax_ids.setdefault(seq_id.split(":")[0], None)
        seq_ids.append(seq_id)
    return carriers


def _connected_groups(links):
    """Connected components (sets) of the items linked by each list in links."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for members in links:
        for m in members[1:]:
            a, b = find(members[0]), find(m)
            if a != b:
                parent[b] = a
    groups = defaultdict(set)
    for x in parent:
        groups[find(x)].add(x)
    return list(groups.values())


def carriers_items(fasta_path, taxonomy):
    """Sequences carried by more than one species with a genus, as (md5, (tax IDs, record IDs))."""
    for md5, (tax_ids, seq_ids) in sequence_carriers(fasta_path).items():
        tax_ids = [t for t in tax_ids if taxonomy.get(t, {}).get("genus")]
        if len(tax_ids) > 1:
            yield md5, (tax_ids, [s for s in seq_ids if s.split(":")[0] in tax_ids])


def identical_reports(fasta_path, taxonomy_path):
    """The rows of identical_16s.tsv and cross_genus_identical.tsv, as lists of dicts."""
    taxonomy = read_taxonomy(taxonomy_path)
    within, cross = [], []
    for md5, (tax_ids, seq_ids) in carriers_items(fasta_path, taxonomy):
        by_genus = defaultdict(list)
        for t in tax_ids:
            by_genus[strip_gtdb_suffixes(taxonomy[t]["genus"])].append(t)
        # Species of one genus that carry the sequence are linked, even when a
        # species of another genus (a contaminant contig) carries it too
        within.extend(members for members in by_genus.values() if len(members) > 1)
        if len(by_genus) > 1:
            cross.append({
                "md5": md5, "species_count": len(tax_ids), "sequence_ids": ";".join(seq_ids),
                "tax_ids": ";".join(tax_ids),
                **{col: ";".join(sorted({taxonomy[t][rank] for t in tax_ids}))
                   for col, rank in (("species", "species"), ("genera", "genus"),
                                     ("families", "family"), ("phyla", "phylum"))},
            })

    shares = defaultdict(set)
    for tax_ids in within:
        for t in tax_ids:
            shares[t].update(u for u in tax_ids if u != t)
    groups = sorted(_connected_groups(within),
                    key=lambda g: (-len(g), min(taxonomy[t]["species"] for t in g)))
    rows = []
    for gi, group in enumerate(groups, 1):
        for t in sorted(group, key=lambda t: (taxonomy[t]["species"], int(t) if t.isdigit() else 0, t)):
            rows.append({
                "tax_id": t, "species": taxonomy[t]["species"], "genus": taxonomy[t]["genus"],
                "group": f"G{gi}", "group_size": len(group),
                "identical_16s_species": "; ".join(sorted({taxonomy[u]["species"] for u in shares[t]})),
            })
    return rows, cross


def write_tsv(rows, columns, path):
    with open(path, "w", newline="", encoding="utf8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_reports(fasta_path, taxonomy_path, groups_path, cross_genus_path=None):
    """Write identical_16s.tsv (and cross_genus_identical.tsv when its path is given); return the rows."""
    rows, cross = identical_reports(fasta_path, taxonomy_path)
    write_tsv(rows, GROUP_COLUMNS, groups_path)
    if cross_genus_path:
        write_tsv(cross, CROSS_GENUS_COLUMNS, cross_genus_path)
    n_groups = len({r["group"] for r in rows})
    print(f"Species sharing an identical 16S within their genus: {len(rows):,} in {n_groups:,} groups ({groups_path})")
    if cross_genus_path:
        print(f"Sequences identical across genera: {len(cross):,} in {sum(c['species_count'] for c in cross):,} "
              f"species entries ({cross_genus_path})")
    return rows, cross


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--fasta", required=True, help="the database's species_taxid.fasta")
    parser.add_argument("--taxonomy", required=True, help="the database's taxonomy.tsv")
    parser.add_argument("--output", required=True, help="output identical_16s.tsv")
    parser.add_argument("--cross-genus", help="also write cross_genus_identical.tsv here")
    args = parser.parse_args()
    write_reports(args.fasta, args.taxonomy, args.output, args.cross_genus)


if __name__ == "__main__":
    main()
