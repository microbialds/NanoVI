#!/usr/bin/env python3
"""Comparator tools (Phase 7a, item 4; decision D8).

  emu          Emu 3.6.2 (quay.io/biocontainers/emu:3.6.2--hdfd78af_0; minimap2 2.30),
               defaults (-N 50, -K 500M, map-ont, --min-abundance 1e-4,
               --max-align-len 2000), --keep-counts. Database 'emu_default' (Emu's
               prebuilt rrnDB v5.6 + NCBI 16S RefSeq, OSF 56uf7 emu.tar) or
               'nanovi_exact' (NanoVI's GTDB r226 database, which has Emu's layout)
  kraken2      Kraken2 2.1.3 + Bracken 3.1 on EPI2ME's ncbi_16s_18s database
               (the one wf-16s downloads; Bracken's 1000-mer distribution, species
               level, threshold 10 reads, as wf-16s)
  wf16s        EPI2ME wf-16s v1.6.1 with --classifier minimap2 or kraken2, its
               default database (ncbi_16s_18s) and settings (minimap2 mode: identity
               >= 95%, reference coverage >= 90%; reads 800-2,000 bp), species rank
               (--taxonomic_rank S; the default is genus), Nextflow 25.04.8 and Docker
  nanoclust    NanoCLUST (genomicsITER/NanoCLUST, commit 9364ddc, DSL1, Nextflow 22.10.8,
               its Docker images), defaults (fastp -q 8, reads 1,400-1,700 bp, UMAP +
               HDBSCAN, min_cluster_size 50, canu/racon/medaka consensus, BLAST against
               NCBI 16S_ribosomal_RNA + taxdb, downloaded 2026-10-23); abundances are the
               clusters' read shares by BLAST top-hit species (<sample>.nanoclust_out.txt,
               so the unipept lookup of its last step is not needed). Real Zymo mocks only (D20)

Names are mapped to GTDB r226 (benchmarks/names.py) before comparison; the
unclassified and unassigned rows are left out of the abundances.
"""
import argparse
import gzip
import os
import subprocess
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, write_tsv
from names import load as load_crosswalk, to_gtdb

EMU_IMAGE = "quay.io/biocontainers/emu:3.6.2--hdfd78af_0"
KRAKEN_IMAGE = "quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0"
BRACKEN_IMAGE = "quay.io/biocontainers/bracken:3.1--h9948957_0"
WF16S = ("epi2me-labs/wf-16s", "v1.6.1")
EMU_DBS = {"emu_default": BENCH / "downloads" / "comparator_db" / "emu_default", "nanovi_exact": BENCH / "db" / "exact"}
KRAKEN_DB = BENCH / "downloads" / "comparator_db" / "ncbi_16s_18s" / "kraken2"
NANOCLUST = BENCH / "downloads" / "nanoclust"
BLAST_IMAGE = "quay.io/biocontainers/blast:2.17.0--h66d330f_0"


def run_emu(reads, db, out_dir, threads=4, N=None):
    """Emu on reads (.fastq.gz) or on a SAM, which Emu takes in place of the reads (it then
    skips minimap2). Each sample gets its own working folder, so runs can go in parallel."""
    out_dir = Path(out_dir)
    name = Path(reads).name.replace(".fastq.gz", "").removesuffix(".sam")
    result = out_dir / f"{name}_rel-abundance.tsv"
    extra = f"--N {N} " if N else ""
    tmp = out_dir / f"tmp_{name}"
    cmd = (f"emu abundance {reads} --db {EMU_DBS[db]} --threads {threads} --keep-counts {extra}"
           f"--output-dir {tmp} --output-basename {name} && mv {tmp}/{name}_rel-abundance.tsv {{OUT}}")
    return harness._run(result, harness.docker(EMU_IMAGE, cmd, [Path(reads).parent, EMU_DBS[db], out_dir]))


def parse_emu(path, ncbi_names):
    t = pd.read_csv(path, sep="\t", dtype={"tax_id": str})
    t = t[~t["tax_id"].isin(["unassigned", "unmapped", "mapped_unclassified"])]
    est = {}
    for sp, ab in zip(t["species"], t["abundance"]):
        name = to_gtdb(sp, ncbi_names) if isinstance(sp, str) else "unnamed"
        est[name] = est.get(name, 0.0) + float(ab)
    total = sum(est.values())
    return {k: v / total for k, v in est.items()}


def run_kraken_bracken(reads, out_dir, threads=4):
    out_dir = Path(out_dir)
    name = Path(reads).name.replace(".fastq.gz", "")
    report = out_dir / f"{name}.kraken2.report"
    harness._run(report, harness.docker(
        KRAKEN_IMAGE, f"kraken2 --db {KRAKEN_DB} --threads {threads} --report {{OUT}} --output /dev/null {reads}",
        [Path(reads).parent, KRAKEN_DB, out_dir]))
    result = out_dir / f"{name}.bracken.tsv"
    return harness._run(result, harness.docker(
        BRACKEN_IMAGE, f"bracken -d {KRAKEN_DB} -i {report} -o {{OUT}} -w {out_dir}/{name}.bracken.report "
                       f"-r 1000 -l S -t 10", [KRAKEN_DB, out_dir]))


def parse_bracken(path, ncbi_names):
    t = pd.read_csv(path, sep="\t")
    est = {}
    for sp, n in zip(t["name"], t["new_est_reads"]):
        name = to_gtdb(sp, ncbi_names)
        est[name] = est.get(name, 0.0) + float(n)
    total = sum(est.values())
    return {k: v / total for k, v in est.items()} if total else {}


def run_wf16s(fastq_dir, out_dir, classifier, threads=8):
    """One wf-16s run over a directory with one sub-directory of reads per sample."""
    out_dir = Path(out_dir)
    done = out_dir / "abundance_table_species.tsv"
    if done.exists():
        return done
    work = out_dir.with_name(out_dir.name + "_work")
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "NXF_VER": os.environ.get("WF16S_NXF_VER", "25.04.8")}
    subprocess.run(["nextflow", "-q", "run", WF16S[0], "-r", WF16S[1], "-profile", "standard",
                    "-w", str(work), "--fastq", str(fastq_dir), "--classifier", classifier,
                    "--database_set", "ncbi_16s_18s", "--taxonomic_rank", "S", "--threads", str(threads),
                    "--store_dir", str(BENCH / "downloads" / "comparator_db" / "wf16s_store"),
                    "--out_dir", str(out_dir), "--include_read_assignments"],
                   check=True, env=env, stdout=open(str(out_dir) + ".log", "w"), stderr=subprocess.STDOUT)
    if not done.exists():
        raise RuntimeError(f"wf-16s wrote no species table in {out_dir}")
    return done


def parse_wf16s(path, sample, ncbi_names):
    t = pd.read_csv(path, sep="\t")
    col = sample if sample in t.columns else next(c for c in t.columns if c.startswith(sample))
    est = {}
    for tax, n in zip(t["tax"], t[col]):
        sp = str(tax).split(";")[-1].strip()
        if not n or sp.lower() in ("unknown", "unclassified", "") or sp.startswith("Unknown"):
            continue
        name = to_gtdb(sp, ncbi_names)
        est[name] = est.get(name, 0.0) + float(n)
    total = sum(est.values())
    return {k: v / total for k, v in est.items()} if total else {}


def run_nanoclust(reads, out_dir, cpus=16):
    """One NanoCLUST run on one sample (gzipped FASTQ; NanoCLUST names it after the file)."""
    out_dir = Path(out_dir)
    name = Path(reads).name.replace(".fastq.gz", "")
    table = out_dir / name / f"{name}.nanoclust_out.txt"
    if table.exists():
        return table
    out_dir.mkdir(parents=True, exist_ok=True)
    fastq = out_dir / f"{name}.fastq"
    # header lines keep the read ID only: its k-mer table splits on tabs, and Dorado's
    # header tags (BC:Z:..., qs:f:...) become extra columns that break the clustering step
    with gzip.open(reads, "rt") as fin, open(fastq, "w") as fout:
        for i, line in enumerate(fin):
            fout.write(line.split(maxsplit=1)[0] + "\n" if i % 4 == 0 else line)
    config = out_dir / "local.config"
    # its BLAST 2.10 cannot read the current NCBI 16S database (GI overflow): a newer BLAST
    # runs the same command; the samples share the output folder, so its reports are rewritten
    config.write_text(f"executor {{ cpus = {cpus} }}\n"
                      f"process {{ withName: consensus_classification {{ container = '{BLAST_IMAGE}' }} }}\n"
                      + "".join(f"{r} {{ overwrite = true }}\n" for r in ("trace", "timeline", "report", "dag")))
    env = {**os.environ, "NXF_VER": "22.10.8"}
    with open(out_dir / f"{name}.log", "w") as log:
        # its last steps look names up online (unipept); the table before them is all we parse
        subprocess.run(["nextflow", "run", str(NANOCLUST / "src" / "main.nf"), "-dsl1", "-profile", "docker",
                        "-c", str(config), "-w", str(out_dir / "work"), "--reads", str(fastq),
                        "--db", str(NANOCLUST / "db" / "16S" / "16S_ribosomal_RNA"),
                        "--tax", str(NANOCLUST / "db" / "taxdb"), "--outdir", str(out_dir)],
                       env=env, stdout=log, stderr=subprocess.STDOUT, cwd=out_dir)
    fastq.unlink()
    if not table.exists():
        raise RuntimeError(f"NanoCLUST wrote no cluster table for {name} (see {out_dir}/{name}.log)")
    return table


def parse_nanoclust(path, ncbi_names):
    t = pd.read_csv(path, sep=";")
    est = {}
    for sp, n in zip(t["sciname"], t["reads_in_cluster"]):
        # BLAST 2.17 quotes the names in -outfmt 10
        name = to_gtdb(" ".join(sp.strip('"').split()[:2]), ncbi_names) if isinstance(sp, str) else "unclassified"
        est[name] = est.get(name, 0.0) + float(n)
    total = sum(est.values())
    return {k: v / total for k, v in est.items()} if total else {}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tool", choices=["emu", "kraken2", "wf16s"], required=True)
    ap.add_argument("--reads", nargs="+", help="FASTQ files (emu, kraken2)")
    ap.add_argument("--db", choices=list(EMU_DBS), default="emu_default")
    ap.add_argument("--fastq-dir", help="wf16s: directory of per-sample sub-directories")
    ap.add_argument("--classifier", choices=["minimap2", "kraken2"], default="minimap2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.tool == "emu":
        for r in args.reads:
            run_emu(r, args.db, out, args.threads)
    elif args.tool == "kraken2":
        for r in args.reads:
            run_kraken_bracken(r, out, args.threads)
    else:
        run_wf16s(args.fastq_dir, out, args.classifier, args.threads)


if __name__ == "__main__":
    main()
