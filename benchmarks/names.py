#!/usr/bin/env python3
"""NCBI -> GTDB species crosswalk (Phase 7a).

The comparators' databases (Emu's default, EPI2ME wf-16s, Kraken2) and the mock
communities' truth tables use NCBI names; NanoVI's database and the simulations use
GTDB r226 names with suffixes stripped (build_gtdb_db.py). Every NCBI species is
mapped to the GTDB species that most of its genomes are placed in, from the GTDB
r226 metadata (ncbi_taxonomy and gtdb_taxonomy of every genome, bacteria and
archaea). The share of genomes behind the majority is kept, so that species
split across several GTDB species can be reported.

Writes <bench>/names/ncbi_to_gtdb.tsv: ncbi_species, gtdb_species (stripped),
genomes, majority_share, n_gtdb_species.
"""
import pandas as pd

from common import BENCH, GTDB, write_tsv
from build_gtdb_db import parse_gtdb_lineage


def crosswalk():
    frames = []
    for f in ("bac120_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv.gz"):
        frames.append(pd.read_csv(GTDB / f, sep="\t", usecols=["accession", "ncbi_taxonomy", "gtdb_taxonomy"],
                                  dtype=str, low_memory=False))
    m = pd.concat(frames)
    m["ncbi_species"] = m["ncbi_taxonomy"].str.extract(r"s__([^;]*)$")[0].str.strip()
    m["gtdb_species"] = m["gtdb_taxonomy"].map(lambda t: parse_gtdb_lineage(t)["species"])
    m = m[(m["ncbi_species"].fillna("") != "") & (m["gtdb_species"] != "")]
    counts = m.groupby(["ncbi_species", "gtdb_species"]).size().rename("n").reset_index()
    totals = counts.groupby("ncbi_species")["n"].agg(["sum", "size"])
    best = counts.sort_values(["ncbi_species", "n", "gtdb_species"], ascending=[True, False, True]) \
        .drop_duplicates("ncbi_species").set_index("ncbi_species")
    out = pd.DataFrame({
        "gtdb_species": best["gtdb_species"],
        "genomes": totals["sum"],
        "majority_share": (best["n"] / totals["sum"]).round(4),
        "n_gtdb_species": totals["size"],
    }).reset_index()
    return out


def load():
    t = pd.read_csv(BENCH / "names" / "ncbi_to_gtdb.tsv", sep="\t")
    return dict(zip(t["ncbi_species"], t["gtdb_species"]))


def to_gtdb(name, table):
    """GTDB name of an NCBI species; names with no genome in GTDB are kept as they are."""
    if name in table:
        return table[name]
    # NCBI names in brackets ('[Clostridium] scindens') appear without brackets in some sources
    alt = name.replace("[", "").replace("]", "")
    return table.get(alt, name)


def main():
    out = crosswalk()
    write_tsv(out, BENCH / "names" / "ncbi_to_gtdb.tsv")
    print(f"{len(out)} NCBI species; {int((out['n_gtdb_species'] > 1).sum())} placed in more than one GTDB species")


if __name__ == "__main__":
    main()
