#!/usr/bin/env python3
"""Database variants for the M3 arm (decision D9).

  exact  build_gtdb_db.py on ssu_all (the database the pipeline uses today:
         exact duplicates removed within a species)
  c995   exact, then each species' sequences clustered at 99.5% identity with
         vsearch --cluster_fast; the centroids (the longest sequence of each
         cluster) are kept. Same taxonomy.tsv and tax IDs as exact.
  reps   build_gtdb_db.py on GTDB's ssu_reps (bac120 + ar53): one 16S per
         species-representative genome. Its own tax IDs (names sorted).

exact and reps are built by bin/build_gtdb_db.py directly (see methods_2026.md);
this script builds c995 from exact, and can remove species from a database for
the clade-exclusion runs (--exclude).
"""
import argparse
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import common  # noqa: F401  (puts bin/ on sys.path)
from build_gtdb_db import fasta_records

VSEARCH = os.environ.get("VSEARCH", "vsearch")


def read_db(db_dir):
    by_taxid = defaultdict(list)
    for header, seq in fasta_records(str(Path(db_dir) / "species_taxid.fasta")):
        by_taxid[int(header[1:].split(":")[0])].append((header[1:], seq))
    return by_taxid


def cluster_species(records, identity):
    """Centroid headers of one species' sequences at the given identity."""
    if len(records) == 1:
        return [records[0][0]]
    with tempfile.TemporaryDirectory() as tmp:
        fin, fout = Path(tmp) / "in.fa", Path(tmp) / "centroids.fa"
        fin.write_text("".join(f">{h}\n{s}\n" for h, s in records))
        subprocess.run([VSEARCH, "--cluster_fast", str(fin), "--id", str(identity), "--qmask", "none",
                        "--threads", "1", "--centroids", str(fout), "--quiet"], check=True)
        return [line[1:].split()[0] for line in fout.read_text().splitlines() if line.startswith(">")]


def write_db(by_taxid, keep_headers, src_dir, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out_dir / "species_taxid.fasta", "w") as fh:
        for taxid in sorted(by_taxid):
            for h, s in by_taxid[taxid]:
                if h in keep_headers:
                    fh.write(f">{h}\n{s}\n")
                    n += 1
    shutil.copy(Path(src_dir) / "taxonomy.tsv", out_dir / "taxonomy.tsv")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", required=True, help="source database directory (exact)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--identity", type=float, default=0.995)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--exclude", help="file of tax IDs to remove instead of clustering (clade exclusion)")
    args = ap.parse_args()

    by_taxid = read_db(args.src)
    if args.exclude:
        drop = {int(x) for x in Path(args.exclude).read_text().split()}
        keep = {h for t, recs in by_taxid.items() if t not in drop for h, _ in recs}
    else:
        taxids = sorted(by_taxid)
        with ProcessPoolExecutor(args.threads) as pool:
            results = pool.map(cluster_species, [by_taxid[t] for t in taxids],
                               [args.identity] * len(taxids), chunksize=64)
            keep = {h for hs in results for h in hs}
    n = write_db(by_taxid, keep, args.src, args.out)
    print(f"{sum(map(len, by_taxid.values()))} sequences in, {n} written to {args.out}")


if __name__ == "__main__":
    main()
