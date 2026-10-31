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
'Pseudomonas'), so genera are grouped the same way whether or not a GTDB
database was built with --merge-gtdb-suffixes. The species differ: with the
suffixes kept, GTDB siblings such as 'Listeria monocytogenes' and 'Listeria
monocytogenes_B' are separate species of one genus and can form groups (GTDB
r226: 5,215 species in 1,513 groups, against 4,335 in 1,264 merged). Species with no genus in taxonomy.tsv are
left out of both files. build_gtdb_db.py and 'main.py build-database' write
both files into the database directory; the pipeline builds identical_16s.tsv
for a database that lacks it (IDENTICAL_16S) and names the sharing species in
the identical_16s_species column of each sample's table.

With a cross-genus rule other than 'keep' (--cross-genus; F56), the builders
also clean the database: the rule removes entries of the sequences in
cross_genus_identical.tsv from species_taxid.fasta, which is rewritten in place,
and identical_16s.tsv is built from what is left. cross_genus_identical.tsv
still lists the sequences as built, since they are what the rule acted on.

  remove       every copy
  best-genus   each sequence kept only in the genus with the most entries
               carrying it (the species of that genus that carry it); a tie
               removes every copy
  most-copies  each sequence kept only in the genus whose entries stand for the
               most copies in the input (a builder deduplicates identical
               sequences of a species into one entry; for GTDB a copy is about
               one genome), the genus with the most database entries breaking a
               tie; a tie on both removes every copy

Two more files list what the rule did:

  cross_genus_removed.tsv    one row per entry removed, with the genus the
                             sequence was kept in ('' when no genus kept it)
  cross_genus_emptied.tsv    the species left with no sequence (they stay in
                             taxonomy.tsv, as species without a long enough
                             sequence do)
"""
import argparse
import csv
import gzip
import hashlib
import os
import re
from collections import defaultdict

# GTDB clade suffix ('Pseudomonas_E', 'Clostridium_AQ'): one or two capitals that end a
# word (GTDB r226 has no longer ones), so other databases' names ('Lachnospiraceae_NK4A136_group',
# 'Family_XIII') stay whole; build_gtdb_db.py strips it with strip_gtdb_suffixes too
GTDB_SUFFIX = re.compile(r"_[A-Z]{1,2}(?=\s|$)")

GROUP_COLUMNS = ["tax_id", "species", "genus", "group", "group_size", "identical_16s_species"]
CROSS_GENUS_COLUMNS = ["md5", "species_count", "sequence_ids", "tax_ids", "species", "genera", "families", "phyla"]
CROSS_GENUS_RULES = ("keep", "remove", "best-genus", "most-copies")
REMOVED_COLUMNS = ["sequence_id", "tax_id", "species", "genus", "md5", "rule", "kept_in"]
EMPTIED_COLUMNS = ["tax_id", "species", "genus", "sequences_removed"]


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


def _unique_best(genera, *scores):
    """The genus with the highest first score, each next score breaking a tie; '' if still tied."""
    best = list(genera)
    for score in scores:
        top = max(score(g) for g in best)
        best = [g for g in best if score(g) == top]
        if len(best) == 1:
            return best[0]
    return ""


def cross_genus_removals(cross, taxonomy, rule, copies=None, genus_entries=None):
    """Entries a cross-genus rule removes: {sequence ID: genus the sequence was kept in, or ''}.

    cross holds the rows of cross_genus_identical.tsv, taxonomy is read_taxonomy's
    dict. Genera are compared with GTDB's suffixes stripped, as in the reports.
    most-copies needs copies (sequence ID -> input copies; missing IDs count 1)
    and genus_entries (genus -> database entries) for its tie-break.
    """
    if rule not in CROSS_GENUS_RULES:
        raise ValueError(f"unknown cross-genus rule {rule!r}; choose one of {', '.join(CROSS_GENUS_RULES)}")
    removed = {}
    if rule == "keep":
        return removed
    copies = copies or {}
    for row in cross:
        by_genus = defaultdict(list)
        for seq_id in row["sequence_ids"].split(";"):
            by_genus[strip_gtdb_suffixes(taxonomy[seq_id.split(":")[0]]["genus"])].append(seq_id)
        kept = ""
        if rule == "best-genus":
            kept = _unique_best(by_genus, lambda g: len(by_genus[g]))
        elif rule == "most-copies":
            kept = _unique_best(by_genus, lambda g: sum(copies.get(s, 1) for s in by_genus[g]),
                                lambda g: genus_entries.get(g, 0))
        for genus, seq_ids in by_genus.items():
            if genus != kept:
                removed.update(dict.fromkeys(seq_ids, kept))
    return removed


def apply_cross_genus_rule(fasta_path, taxonomy_path, cross, rule, out_dir, copies=None):
    """Remove the entries a rule picks from fasta_path (rewritten in place) and write its two reports.

    Writes cross_genus_removed.tsv and cross_genus_emptied.tsv into out_dir and
    returns (entries removed, species left with no sequence).
    """
    taxonomy = read_taxonomy(taxonomy_path)
    entries = defaultdict(int)
    for header, _ in fasta_records(fasta_path):
        entries[header[1:].split()[0].split(":")[0]] += 1
    genus_entries = defaultdict(int)
    for t, n in entries.items():
        genus_entries[strip_gtdb_suffixes(taxonomy.get(t, {}).get("genus", ""))] += n
    removed = cross_genus_removals(cross, taxonomy, rule, copies, genus_entries)
    md5_of = {s: row["md5"] for row in cross for s in row["sequence_ids"].split(";")}

    tmp = fasta_path + ".tmp"
    lost = defaultdict(int)
    with open(tmp, "w", encoding="utf8") as fh:
        for header, seq in fasta_records(fasta_path):
            seq_id = header[1:].split()[0]
            if seq_id in removed:
                lost[seq_id.split(":")[0]] += 1
            else:
                fh.write(f"{header}\n{seq}\n")
    os.replace(tmp, fasta_path)

    write_tsv([{"sequence_id": s, "tax_id": s.split(":")[0], "species": taxonomy[s.split(":")[0]]["species"],
                "genus": taxonomy[s.split(":")[0]]["genus"], "md5": md5_of[s], "rule": rule, "kept_in": kept}
               for s, kept in removed.items()], REMOVED_COLUMNS, os.path.join(out_dir, "cross_genus_removed.tsv"))
    emptied = sorted((t for t, n in lost.items() if n == entries[t]),
                     key=lambda t: (taxonomy[t]["species"], int(t) if t.isdigit() else 0, t))
    write_tsv([{"tax_id": t, "species": taxonomy[t]["species"], "genus": taxonomy[t]["genus"],
                "sequences_removed": lost[t]} for t in emptied],
              EMPTIED_COLUMNS, os.path.join(out_dir, "cross_genus_emptied.tsv"))
    print(f"Cross-genus rule '{rule}': {len(removed):,} entries of {len(lost):,} species removed; "
          f"{len(emptied):,} species left with no sequence (cross_genus_removed.tsv, cross_genus_emptied.tsv)")
    return len(removed), len(emptied)


def write_tsv(rows, columns, path):
    with open(path, "w", newline="", encoding="utf8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_reports(fasta_path, taxonomy_path, groups_path, cross_genus_path=None, cross_genus_rule="keep",
                  copies=None):
    """Write identical_16s.tsv (and cross_genus_identical.tsv when its path is given); return the rows.

    A cross_genus_rule other than 'keep' first cleans fasta_path (apply_cross_genus_rule,
    reports next to groups_path; copies as there), so identical_16s.tsv describes the
    cleaned database, while cross_genus_identical.tsv lists the sequences as built.
    """
    rows, cross = identical_reports(fasta_path, taxonomy_path)
    if cross_genus_path:
        write_tsv(cross, CROSS_GENUS_COLUMNS, cross_genus_path)
    if cross_genus_rule != "keep":
        apply_cross_genus_rule(fasta_path, taxonomy_path, cross, cross_genus_rule,
                               os.path.dirname(os.path.abspath(groups_path)), copies)
        rows, _ = identical_reports(fasta_path, taxonomy_path)
    write_tsv(rows, GROUP_COLUMNS, groups_path)
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
