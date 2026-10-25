#!/usr/bin/env python3
"""Seeded simulated communities with known truth (M9).

One community of 20 species: five near-neighbour pairs with both members present
and ten unrelated species. Each species is represented by one GTDB genome (chosen
once, seeded), and reads are simulated from all of that genome's 16S copies of at
least 1,400 bp in ssu_all, so intragenomic variation is included. The genome's
copies are in the exact database (the database has every distinct sequence of a
species), not necessarily in c995 (a centroid within 0.5%) or reps (one genome).

Settings: depth {1k, 10k, 50k} reads x chemistry {r941: Badread nanopore2020 models,
identity 94,98,3 (~6% errors); r1041: nanopore2023 models, 99,100,0.5 (~1%)} x shape
{even: 5% each; log: abundances 10^(-4 i / 19) normalised, spanning four orders of
magnitude, pair members ranked apart} x seeds. Per sample, species counts are a
multinomial draw from the design abundances (seeded), then each species' reads come
from one Badread call (whole references: fragment length 20 kb on 16S copies;
adapters and glitches at Badread's defaults; no junk, random or chimeric reads).

Writes <bench>/sims/<sample>/reads.fastq.gz, truth.tsv (species, design abundance,
reads) and read_truth.tsv.gz (read, species, genome copy).
"""
import argparse
import gzip
import os
import random
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from common import BENCH, SSU_ALL, write_tsv
from build_gtdb_db import fasta_records, parse_ssu_header

BADREAD = os.environ.get("BADREAD", "badread")

PAIRS = [
    ("Escherichia coli", "Escherichia fergusonii"),
    ("Klebsiella pneumoniae", "Klebsiella variicola"),
    ("Streptococcus pneumoniae", "Streptococcus mitis"),
    ("Bacillus subtilis", "Bacillus spizizenii"),
    ("Staphylococcus aureus", "Staphylococcus argenteus"),
]
OTHERS = [
    "Pseudomonas aeruginosa", "Enterococcus faecalis", "Listeria monocytogenes", "Salmonella enterica",
    "Limosilactobacillus fermentum", "Bacteroides fragilis", "Akkermansia muciniphila",
    "Clostridioides difficile", "Bifidobacterium longum", "Acinetobacter baumannii",
]
# Rank order for the log-staggered shape: the second member of each pair sits 5-9 ranks
# below the first (about 10-50 times less abundant)
LOG_ORDER = [
    "Escherichia coli", "Pseudomonas aeruginosa", "Staphylococcus aureus", "Klebsiella pneumoniae",
    "Enterococcus faecalis", "Streptococcus pneumoniae", "Bacillus subtilis", "Escherichia fergusonii",
    "Listeria monocytogenes", "Salmonella enterica", "Klebsiella variicola", "Bacteroides fragilis",
    "Staphylococcus argenteus", "Limosilactobacillus fermentum", "Streptococcus mitis",
    "Akkermansia muciniphila", "Bacillus spizizenii", "Clostridioides difficile",
    "Bifidobacterium longum", "Acinetobacter baumannii",
]
CHEMISTRY = {
    "r941": ["--error_model", "nanopore2020", "--qscore_model", "nanopore2020", "--identity", "94,98,3"],
    "r1041": ["--error_model", "nanopore2023", "--qscore_model", "nanopore2023", "--identity", "99,100,0.5"],
}
COMMON_ARGS = ["--length", "20000,10", "--junk_reads", "0", "--random_reads", "0", "--chimeras", "0"]
MIN_COPY_LENGTH = 1400


def community_species():
    return [s for pair in PAIRS for s in pair] + OTHERS


def design_abundances(shape):
    species = community_species()
    if shape == "even":
        return {s: 1 / len(species) for s in species}
    w = 10 ** (-4 * np.arange(len(LOG_ORDER)) / (len(LOG_ORDER) - 1))
    w /= w.sum()
    return dict(zip(LOG_ORDER, w))


def choose_genomes(seed=7):
    """One genome per species: RefSeq first, then most 16S copies >= 1,400 bp, then seeded choice.

    Copies whose sequence also occurs in another genus (a contaminant contig in some
    assembly; db_audit.py) disqualify a genome.
    """
    audit = pd.read_csv(BENCH / "audit" / "ssu_all" / "sequences.tsv.gz", sep="\t")
    genus_spread = audit.drop_duplicates(["species", "md5"]).groupby("md5")["genus"].nunique()
    bad_md5 = set(genus_spread[genus_spread > 1].index)
    rng = random.Random(seed)
    chosen = {}
    for sp in community_species():
        rows = audit[(audit["species"] == sp) & (audit["length"] >= MIN_COPY_LENGTH)]
        bad_genomes = set(rows.loc[rows["md5"].isin(bad_md5), "genome"])
        g = rows[~rows["genome"].isin(bad_genomes)].groupby("genome").size().rename("copies").reset_index()
        g["refseq"] = g["genome"].str.startswith("RS_")
        best = g[g["refseq"] == g["refseq"].max()]
        best = best[best["copies"] == best["copies"].max()]
        chosen[sp] = rng.choice(sorted(best["genome"]))
    return chosen


def write_references(genomes, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    by_genome = {g: sp for sp, g in genomes.items()}
    copies = {sp: [] for sp in genomes}
    for header, seq in fasta_records(str(SSU_ALL)):
        genome, lineage = parse_ssu_header(header)
        if genome in by_genome and len(seq) >= MIN_COPY_LENGTH:
            copies[by_genome[genome]].append((header[1:].split()[0], seq))
    rows = []
    for i, sp in enumerate(community_species()):
        path = out_dir / f"sp{i:02d}.fasta"
        with open(path, "w") as fh:
            for k, (sid, seq) in enumerate(copies[sp]):
                fh.write(f">sp{i:02d}_c{k}\n{seq}\n")
                rows.append((f"sp{i:02d}", sp, genomes[sp], f"sp{i:02d}_c{k}", sid, len(seq)))
    df = pd.DataFrame(rows, columns=["sp", "species", "genome", "copy", "gtdb_ssu_id", "length"])
    write_tsv(df, out_dir / "copies.tsv")
    return df


def run_badread(args):
    """One Badread call; its output is kept, so an interrupted run resumes where it stopped."""
    ref, n_reads, chem, seed, out = args
    if Path(out).exists():
        return out
    mean_len = np.mean([len(s) for _, s in fasta_records(str(ref))])
    quantity = int(n_reads * mean_len * 1.25) + 5000
    tmp = Path(str(out) + ".tmp")
    with open(tmp, "w") as fh, open(str(out) + ".log", "w") as log:
        subprocess.run([BADREAD, "simulate", "--reference", str(ref), "--quantity", str(quantity),
                        "--seed", str(seed), *CHEMISTRY[chem], *COMMON_ARGS],
                       stdout=fh, stderr=log, check=True)
    tmp.rename(out)
    return out


def fastq_records(path):
    with open(path) as fh:
        while True:
            lines = [fh.readline() for _ in range(4)]
            if not lines[0]:
                return
            yield lines


def sample_name(depth, chem, shape, seed):
    return f"{shape}_{chem}_{depth // 1000}k_s{seed:02d}"


def plan_sample(depth, chem, shape, seed):
    design = design_abundances(shape)
    rng = np.random.default_rng([seed, depth, list(CHEMISTRY).index(chem), ["even", "log"].index(shape)])
    species = community_species()
    counts = rng.multinomial(depth, [design[s] for s in species])
    return design, dict(zip(species, counts.tolist()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--depths", default="1000,10000,50000")
    ap.add_argument("--chemistries", default="r941,r1041")
    ap.add_argument("--shapes", default="even,log")
    ap.add_argument("--seeds", default="1-10")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()
    lo, _, hi = args.seeds.partition("-")
    seeds = range(int(lo), int(hi or lo) + 1)

    sims = BENCH / "sims"
    ref_dir = sims / "references"
    if not (ref_dir / "copies.tsv").exists():
        genomes = choose_genomes()
        write_references(genomes, ref_dir)
    copies = pd.read_csv(ref_dir / "copies.tsv", sep="\t")
    sp_of = dict(zip(copies["sp"], copies["species"]))
    species = community_species()

    samples, jobs = [], []
    for depth in map(int, args.depths.split(",")):
        for chem in args.chemistries.split(","):
            for shape in args.shapes.split(","):
                for seed in seeds:
                    name = sample_name(depth, chem, shape, seed)
                    out = sims / name
                    if (out / "reads.fastq.gz").exists():
                        continue
                    (out / "parts").mkdir(parents=True, exist_ok=True)
                    design, counts = plan_sample(depth, chem, shape, seed)
                    samples.append((name, design, counts))
                    for i, sp in enumerate(species):
                        if counts[sp]:
                            jobs.append((ref_dir / f"sp{i:02d}.fasta", counts[sp], chem,
                                         seed * 1000 + i, out / "parts" / f"sp{i:02d}.fastq"))
    print(f"{len(samples)} samples, {len(jobs)} Badread calls")
    with ProcessPoolExecutor(args.threads) as pool:
        for i, _ in enumerate(pool.map(run_badread, jobs, chunksize=1), 1):
            if i % 100 == 0:
                print(f"  {i}/{len(jobs)} Badread calls done", flush=True)

    for name, design, counts in samples:
        out = sims / name
        records, truth = [], []
        rng = random.Random(name)
        for i, sp in enumerate(species):
            if not counts[sp]:
                continue
            part = list(fastq_records(out / "parts" / f"sp{i:02d}.fastq"))
            if len(part) < counts[sp]:
                raise RuntimeError(f"{name} {sp}: Badread gave {len(part)} reads, {counts[sp]} needed")
            for rec in part[:counts[sp]]:
                read_id, source = rec[0][1:].split()[:2]
                truth.append((read_id, sp, source.split(",")[0]))
                records.append(rec)
        rng.shuffle(records)
        with gzip.open(out / "reads.fastq.gz", "wt", compresslevel=6) as fh:
            for rec in records:
                fh.writelines(rec)
        write_tsv(pd.DataFrame(truth, columns=["read", "species", "copy"]), out / "read_truth.tsv.gz")
        write_tsv(pd.DataFrame({"species": species, "design_abundance": [design[s] for s in species],
                                "reads": [counts[s] for s in species]}), out / "truth.tsv")
        for p in (out / "parts").iterdir():
            p.unlink()  # parts and their Badread logs
        (out / "parts").rmdir()
    print("done")


if __name__ == "__main__":
    main()
