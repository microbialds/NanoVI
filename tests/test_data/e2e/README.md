# End-to-end test fixture

A small real database and simulated ONT-like reads with known proportions. The
`e2e` nf-test (`tests/main.nf.test`) and the `test` profile (`conf/test.config`)
run the full pipeline on it: fastplong, minimap2 and every Python step.

```
nf-test test tests/main.nf.test --tag e2e          # needs Docker
nextflow run main.nf -profile test,docker          # same data, outputs in test_results/
```

## Contents

| File | What it is |
|---|---|
| `db/species_taxid.fasta` | 11 16S rRNA sequences, one per species, headers `taxid:e2e:n` |
| `db/taxonomy.tsv` | Lineages in the format `bin/build_gtdb_db.py` writes |
| `sources.tsv` | NCBI accession, species and the lineage used for each sequence |
| `reads/e2e_mock.fastq.gz` | 630 simulated reads: 600 16S reads and 30 non-16S reads |
| `expected.tsv` | True proportions and read counts per species |
| `make_reads.py` | The read simulator (fixed seed; the same reads on every Python 3) |

## Database

The sequences are NCBI RefSeq 16S records (BioProject PRJNA33175, mostly type
strains) plus the classic *E. coli* rrnB 16S sequence (J01859.1). Lineages follow
NCBI taxonomy names, written in GTDB header style so that `bin/build_gtdb_db.py`
builds the database exactly as it would from a GTDB SSU file.

| Accession | Species | Length (bp) | In the sample |
|---|---|---|---|
| J01859.1 | *Escherichia coli* | 1,541 | 25 % |
| NR_026331.1 | *Shigella flexneri* ATCC 29903 | 1,488 | 5 % |
| NR_074910.1 | *Salmonella enterica* LT2 | 1,544 | absent |
| NR_117678.1 | *Pseudomonas aeruginosa* DSM 50071 | 1,527 | 10 % |
| NR_028665.1 | *Streptococcus pneumoniae* ATCC 33400 | 1,515 | 15 % |
| NR_028664.1 | *Streptococcus mitis* NS51 | 1,520 | 5 % |
| NR_114782.1 | *Enterococcus faecalis* LMG 7937 | 1,556 | absent |
| NR_104927.1 | *Limosilactobacillus fermentum* CIP 102980 | 1,502 | absent |
| NR_118997.2 | *Staphylococcus aureus* ATCC 12600 | 1,552 | 40 % |
| NR_102783.2 | *Bacillus subtilis* 168 | 1,550 | absent |
| NR_044823.1 | *Listeria monocytogenes* NCTC 10357 | 1,469 | absent |

Two near-neighbour pairs: *E. coli* / *S. flexneri* (99.8 % identity over the
alignment) and *S. pneumoniae* / *S. mitis* (99.3 %). *S. enterica* (97.4 % to
*E. coli*) is a close decoy that is absent from the sample.

To rebuild the database:

```bash
ids=$(tail -n +2 sources.tsv | cut -f1 | paste -sd,)
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ids}&rettype=fasta&retmode=text" > ncbi_16s.fasta
# Rewrite each header as ">ACC~ACC <lineage from sources.tsv>", then:
python3 ../../../bin/build_gtdb_db.py --ssu gtdb_style_ssu.fna --db-name e2e --output-dir db --min-length 1000
```

## Reads

`make_reads.py` draws 600 reads from the database in the proportions
40/25/15/10/5/5 % (exact counts 240/150/90/60/30/30, see `expected.tsv`): a
random 1,200-1,600 bp window of the reference (capped at the reference length),
a random strand, IUPAC codes resolved to a random base, and a per-read error
rate drawn uniformly from 5-8 % (40 % substitutions, 25 % insertions, 35 %
deletions). After indels the read lengths range from 1,178 to 1,599 bp. It then
adds 30 reads of uniformly random sequence that align to nothing, and shuffles
all reads. Quality strings are a constant Phred 12; the pipeline disables
fastplong's quality filter.

```bash
python3 tests/test_data/e2e/make_reads.py
```

## Expected results

The `e2e` test checks that the assigned abundances sum to 1, that every species
at 10 % or more is reported within 5 percentage points of `expected.tsv`, that
every species in the sample is reported (so a near-neighbour pair is not merged
into one species), that no absent species reaches 1 %, and that
`consolidated_output/` holds the four tables. The Unassigned row (30 non-16S
reads) is not checked until F01 is fixed.

Reference output (`-profile test,docker`, which uses the pipeline defaults k = 21 and N = 3; NanoVI 1.0.0 code, October 2026):

| Species | Expected | Reported | estimated_counts |
|---|---|---|---|
| *Staphylococcus aureus* | 0.40 | 0.3996 | 239 |
| *Escherichia coli* | 0.25 | 0.2466 | 147 |
| *Streptococcus pneumoniae* | 0.15 | 0.1499 | 89 |
| *Pseudomonas aeruginosa* | 0.10 | 0.1001 | 60 |
| *Shigella flexneri* | 0.05 | 0.0535 | 32 |
| *Streptococcus mitis* | 0.05 | 0.0503 | 30 |

600 reads assigned, 30 unassigned. No absent species is reported. The table has
no Unassigned row (F01) and the counts sum to 597, not 600 (F08).
