#!/usr/bin/env python3
"""Counts of the two GTDB r226 builds of bin/build_gtdb_db.py (Phase 9, for the Methods).

  merged    --merge-gtdb-suffixes, the database of the Phase 7a benchmark
  suffixes  the default build, GTDB names with their clade suffixes (Phase 7b, F16)

Reads each build's species_taxid.fasta, taxonomy.tsv, identical_16s.tsv and
cross_genus_identical.tsv and writes docs/benchmarks/data/db_builds_r226.tsv.

    python benchmarks/db_builds.py --merged <dir> --suffixes <dir>
"""
import argparse
from pathlib import Path

import pandas as pd

from common import DOCS_DATA, write_tsv


def counts(db):
    db = Path(db)
    seqs, taxa = 0, set()
    with open(db / "species_taxid.fasta") as f:
        for line in f:
            if line.startswith(">"):
                seqs += 1
                taxa.add(line[1:].split(":", 1)[0])
    ident = pd.read_csv(db / "identical_16s.tsv", sep="\t")
    cross = pd.read_csv(db / "cross_genus_identical.tsv", sep="\t", dtype=str)
    carriers = {t for ids in cross["tax_ids"] for t in ids.split(";")}
    return {
        "database sequences": seqs,
        "species with a sequence": len(taxa),
        "species in taxonomy.tsv": len(pd.read_csv(db / "taxonomy.tsv", sep="\t")),
        "species sharing an identical 16S within their genus (identical_16s.tsv)": len(ident),
        "identical-16S groups": ident["group"].nunique(),
        "sequences found in more than one genus (cross_genus_identical.tsv)": len(cross),
        "... of which in more than one family": int((cross["families"].str.count(";") > 0).sum()),
        "species carrying a cross-genus sequence": len(carriers),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--suffixes", required=True)
    args = ap.parse_args()
    m, s = counts(args.merged), counts(args.suffixes)
    df = pd.DataFrame({"metric": list(m), "merged": list(m.values()), "suffixes_kept": [s[k] for k in m]})
    print(df.to_string(index=False))
    write_tsv(df, DOCS_DATA / "db_builds_r226.tsv")


if __name__ == "__main__":
    main()
