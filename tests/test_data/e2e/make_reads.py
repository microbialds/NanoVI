#!/usr/bin/env python3
"""
make_reads.py — simulate the ONT-like 16S reads for the NanoVI e2e fixture.

Draws reads from db/species_taxid.fasta in fixed proportions, adds
substitutions, insertions and deletions at a per-read error rate of 5-8 %,
and appends reads of random sequence that should not align to any 16S
reference. Writes reads/e2e_mock.fastq.gz and expected.tsv.

Only random.Random(SEED).random() is used, so the output is the same on every
Python 3 version, and the gzip header carries no timestamp, so the committed
file is byte-for-byte reproducible:

    python3 tests/test_data/e2e/make_reads.py
"""

import csv
import gzip
import os
import random

SEED = 20261021
HERE = os.path.dirname(os.path.abspath(__file__))

# Reads per species (600 16S reads in total: 40/25/15/10/5/5 %).
# The other five species in the database are absent from the sample.
COMPOSITION = {
    "Staphylococcus aureus": 240,
    "Escherichia coli": 150,
    "Streptococcus pneumoniae": 90,
    "Pseudomonas aeruginosa": 60,
    "Shigella flexneri": 30,
    "Streptococcus mitis": 30,
}
N_NON_16S = 30

MIN_LEN, MAX_LEN = 1200, 1600
MIN_ERR, MAX_ERR = 0.05, 0.08
# Share of errors by type: substitution, insertion, deletion.
ERR_SPLIT = (0.40, 0.25, 0.35)
QUAL_CHAR = "-"  # Phred 12

IUPAC = {
    "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "S": "CG", "W": "AT",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def pick(rng, items):
    return items[int(rng.random() * len(items))]


def uniform_int(rng, lo, hi):
    """Integer in [lo, hi]."""
    return lo + int(rng.random() * (hi - lo + 1))


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif line:
                records[name].append(line.upper())
    return {k: "".join(v) for k, v in records.items()}


def read_taxonomy(path):
    with open(path) as fh:
        return {row["species"]: row["tax_id"] for row in csv.DictReader(fh, delimiter="\t")}


def add_errors(rng, seq, err_rate):
    p_sub = err_rate * ERR_SPLIT[0]
    p_ins = err_rate * ERR_SPLIT[1]
    p_del = err_rate * ERR_SPLIT[2]
    out = []
    for base in seq:
        u = rng.random()
        if u < p_sub:
            out.append(pick(rng, [b for b in "ACGT" if b != base]))
        elif u < p_sub + p_ins:
            out.append(base)
            out.append(pick(rng, "ACGT"))
        elif u < p_sub + p_ins + p_del:
            continue
        else:
            out.append(base)
    return "".join(out)


def simulate_16s_read(rng, ref):
    ref = "".join(pick(rng, IUPAC[b]) if b in IUPAC else b for b in ref)
    length = uniform_int(rng, MIN_LEN, min(MAX_LEN, len(ref)))
    start = uniform_int(rng, 0, len(ref) - length)
    frag = ref[start:start + length]
    if rng.random() < 0.5:
        frag = frag.translate(COMPLEMENT)[::-1]
    err_rate = MIN_ERR + rng.random() * (MAX_ERR - MIN_ERR)
    return add_errors(rng, frag, err_rate)


def random_read(rng):
    return "".join(pick(rng, "ACGT") for _ in range(uniform_int(rng, MIN_LEN, MAX_LEN)))


def main():
    rng = random.Random(SEED)
    refs = read_fasta(os.path.join(HERE, "db", "species_taxid.fasta"))
    taxids = read_taxonomy(os.path.join(HERE, "db", "taxonomy.tsv"))
    ref_by_taxid = {name.split(":")[0]: seq for name, seq in refs.items()}

    reads = []
    for species, n in COMPOSITION.items():
        tid = taxids[species]
        for i in range(n):
            reads.append((f"{tid}_{species.replace(' ', '_')}_{i}", simulate_16s_read(rng, ref_by_taxid[tid])))
    for i in range(N_NON_16S):
        reads.append((f"non16S_{i}", random_read(rng)))

    # Fisher-Yates shuffle with rng.random() only
    for i in range(len(reads) - 1, 0, -1):
        j = int(rng.random() * (i + 1))
        reads[i], reads[j] = reads[j], reads[i]

    os.makedirs(os.path.join(HERE, "reads"), exist_ok=True)
    fastq_path = os.path.join(HERE, "reads", "e2e_mock.fastq.gz")
    with open(fastq_path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            for name, seq in reads:
                gz.write(f"@{name}\n{seq}\n+\n{QUAL_CHAR * len(seq)}\n".encode())

    n_16s = sum(COMPOSITION.values())
    with open(os.path.join(HERE, "expected.tsv"), "w") as fh:
        fh.write("tax_id\tspecies\texpected_abundance\tsimulated_reads\n")
        for species, tid in sorted(taxids.items(), key=lambda kv: int(kv[1])):
            n = COMPOSITION.get(species, 0)
            fh.write(f"{tid}\t{species}\t{n / n_16s:.4f}\t{n}\n")
        fh.write(f"NA\tnon-16S\tNA\t{N_NON_16S}\n")

    print(f"Wrote {len(reads)} reads to {fastq_path}")


if __name__ == "__main__":
    main()
