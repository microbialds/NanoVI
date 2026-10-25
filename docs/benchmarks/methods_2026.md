# NanoVI methods benchmark (Phase 7a, October 2026)

This document records how the 2026 benchmark was run and what it found. It covers the
database audit, simulations with known truth, the ablation of the pipeline's settings,
real mock communities, read filtering, credible intervals, the unknown-taxon prototype,
the comparator tools and runtime. No default changed in this phase; the decisions it
informs are taken in Phase 7b.

Every table cited here is in `docs/benchmarks/data/` and every figure in
`docs/benchmarks/figures/`, with the figure's source data under the same name in
`data/`. The code is in `benchmarks/` (one script per stage, described below); data,
databases, alignments and per-sample outputs are under
`~/pipelines/NanoVI_benchmarks/phase7a` (decision D11) and are not in the repository.

## 1. How the benchmark runs

**Harness.** The pipeline's four abundance steps run outside Nextflow, in the same
containers and with the same commands as the modules (`benchmarks/harness.py`):
FILTER_READS (fastplong 0.2.2: length 500-2,000, quality filter off), ALIGN_READS
(minimap2 2.24 `-ax map-ont -N 3 -p 0.9 -K 500M -k 21` piped into samtools 1.18),
LIKELIHOOD (`bin/read_likelihood.py`) and RUN_VI (`bin/estimate_abundance.py`). Each
step's output is cached under its parameters, so layers that only change the
likelihood or the estimator reuse one alignment. The harness's abundance JSON was
checked identical to `nextflow run main.nf` on the same reads.

**Reference setting.** Exact database, `-N 3`, primary-alignment error profile, Emu's
read score, VB, output floor 1e-3, pruning on: the pipeline's defaults at dev `d8e387b`.

**Compute budget (D20).** At most 16 cores of a shared 80-core host (load 80-250 from
other users during the runs). Aligning against the exact database costs 0.8 CPU-s per
read (Section 3), which made the full D9 grid (>1,500 CPU-h) unaffordable. The design
was cut as follows, and every result below says which part it comes from:

- *One factor at a time (ofat)*: every simulated sample at the reference alignment,
  with one layer changed at a time (error profile, floor, estimator, pruning, score):
  92 samples (1k and 10k reads, 10 seeds; 50k reads, seeds 1-3), 8 arms each.
- *Cell grid*: database × -N × score × profile × estimator × floor on the four
  10k-read seed-1 samples (both chemistries and shapes), 288 arms each.
- Real data subsampled to 10,000 reads; hac basecalls only for the reference and
  quality-filter arms; Emu on the simulations only for the four seed-1 samples;
  NanoCLUST on the real Zymo mocks with an 8-hour budget.

**Environments.** Tools run from BioContainers images with `--network none`;
host-side Python (the pipeline's `bin/` scripts and the analysis) runs in a conda
environment with the pipeline's pins (Python 3.11, numpy 2.4.6, scipy 1.17.1,
pandas 3.0.6, pysam 0.24.1, numba 0.68.0). Lock files are kept beside the data
(`env_*.lock.txt`).

| Tool | Version | Source |
|---|---|---|
| fastplong | 0.2.2 (pipeline), 0.7.0 (D16 arm) | quay.io/biocontainers |
| minimap2 | 2.24 with samtools 1.18 (pipeline), 2.31 (D16 arm) | quay.io/biocontainers |
| Badread | 0.4.2 | bioconda |
| vsearch | 2.32.0 | bioconda |
| seqkit | 2.14.0 | bioconda |
| Emu | 3.6.2 (minimap2 2.30) | quay.io/biocontainers/emu:3.6.2--hdfd78af_0 |
| Kraken2 / Bracken | 2.1.3 / 3.1 | quay.io/biocontainers |
| EPI2ME wf-16s | v1.6.1, Nextflow 25.04.8 | epi2me-labs/wf-16s |
| NanoCLUST | commit 9364ddc, Nextflow 22.10.8 (DSL1) | genomicsITER/NanoCLUST |

**Stages and scripts** (run in this order by `run_phase.sh`, `run_phase2.sh` and
`run_phase3.sh` beside the data):

| Stage | Script | Output |
|---|---|---|
| Database audit and variants | `db_audit.py`, `build_variants.py` | `db_audit_*`, `db_variants_*` |
| Simulations | `simulate.py` | `sims/<sample>/` |
| ofat and cell grids | `run_grid.py --design ofat / cell` | `grid_ofat.tsv`, `grid_cell.tsv` |
| Clade exclusion | `clade_exclusion.py` | `clade_exclusion.tsv` |
| Real mocks | `real.py` | `real_nanovi.tsv` |
| Quality filter on simulations | `quality_sims.py` | `quality_sims.tsv` |
| Credible intervals | `intervals.py` | `intervals.tsv` |
| Comparators | `run_comparators.py` | `comparators.tsv` |
| Runtime | `timing.py` | `timing.tsv` |
| Tables and figures | `analyse.py <analysis>` | `docs/benchmarks/data`, `figures` |

**Metrics** (`benchmarks/evaluate.py`). Species and genus precision, recall and F1;
L1 distance (0-2) and Bray-Curtis dissimilarity to the truth; false-positive count and
abundance. *Group* metrics merge species whose 16S genes are identical (the identical
groups of the audit), the finest level the gene can resolve. On real mocks, L1 is
computed over the species with a known share (the estimate renormalised over them),
recall counts the species expected to give at least 10 reads, and the spike-ins of the
2025.09 release (presence only) count for recall and precision.

## 2. Database audit (M3, M6, W9)

`bin/build_gtdb_db.py` on GTDB r226 `ssu_all` (1,053,664 sequences, 792,620 of at least
900 bp) gives 232,447 sequences for 59,037 species; `taxonomy.tsv` lists 81,173 species
(22,136 without a sequence). The distribution is very uneven
(`db_audit_ssu_all_summary.tsv`, `db_audit_ssu_all_top20.tsv`,
`db_audit_ssu_all_seqs_per_species_hist.tsv`): median 1 sequence per species, 95th
percentile 8, 99th 36, maximum 14,661 (*Escherichia coli*); 39,384 species have one
sequence and the top 1% of species hold 45.7% of the sequences.

Identical 16S across species (`db_audit_ssu_all_identical_groups_top30.tsv`): 3,925
species (6.6%) share an identical sequence with another species of their genus, in
1,238 groups (the largest, 128 *Pseudomonas* species). 1,360 sequences are identical
across genera (1,288 across families): these are almost certainly contaminant contigs in
GTDB assemblies (example: GTDB's only 16S for *Nonomuraea* sp040654455 is a 4,394 bp
contig of RS_GCF_040654455.1 that is 100% identical over 1,550 bp to an *S. aureus* 16S).

Suffix stripping (`db_audit_ssu_all_suffix_merges_top30.tsv`) merges 4,787 GTDB species
into 1,529 names (*Buchnera aphidicola* 111, *Streptococcus mitis* 90, *S. oralis* 75);
26 merged names join species from different GTDB genera. 60,722 names are placeholders
(`sp` plus digits).

**Variants (D9; `db_variants_summary.tsv`).**

| Variant | How | Sequences | Per species: median / p95 / max |
|---|---|---|---|
| exact | build_gtdb_db.py on ssu_all (exact duplicates removed within a species) | 232,447 | 1 / 8 / 14,661 |
| c995 | exact, each species clustered at 99.5% with vsearch `--cluster_fast`, centroids kept | 98,444 | 1 / 3 / 2,615 |
| reps | build_gtdb_db.py on GTDB ssu_reps (bac120 + ar53) | 56,861 | 1 / 1 / 111 |

## 3. Simulated communities (M9)

`simulate.py`: one community of 20 species, five near-neighbour pairs with both members
present (*E. coli* / *E. fergusonii*, *K. pneumoniae* / *K. variicola*,
*S. pneumoniae* / *S. mitis*, *B. subtilis* / *B. spizizenii*, *S. aureus* /
*S. argenteus*) and ten unrelated species. Each species is one GTDB genome; reads come
from all of its 16S copies of at least 1,400 bp, so intragenomic variation is included.
Badread 0.4.2: R9.4.1-like (nanopore2020 error and quality models, identity 94,98,3;
mean read identity 93.9%) and R10.4.1-like (nanopore2023, 99,100,0.5; 99.0%). Even (5%
each) and log-staggered (four orders of magnitude, pair members ranked apart)
communities; 1k, 10k and 50k reads; 10 seeds each, 120 samples (the grids use all 10
at 1k and 10k and seeds 1-3 at 50k, 92 samples). Species counts are a multinomial draw
from the design. Within a chemistry, each species' reads come from one Badread stream
per seed (seed × 1000 + species index) whatever the depth or shape, so a species' reads
at 1k are a prefix of its reads at 10k and 50k with the same seed: replicates within a
condition are independent, but conditions that share a seed share reads.

Two pairs are not identifiable from 16S in GTDB r226: the chosen *S. argenteus* and
*E. fergusonii* genomes carry 16S copies identical to *S. aureus* and *E. coli* copies
(identical groups G72 and G127), so their reads score the same against the partner.
GTDB places *Shigella* in *E. coli*, so the e2e fixture's *E. coli* / *S. flexneri*
pair does not exist in this database.

## 4. Ablation (M1, M3, M4, M7; F18 -N)

### 4.1 Reference setting

On the ofat samples (`ablation_ofat_summary.tsv`), mean over seeds:

| Reads | Chemistry | Shape | Species L1 | Group L1 | Genus L1 | Recall | Precision | FP species | FP abundance |
|---|---|---|---|---|---|---|---|---|---|
| 1k | R10.4.1 | even | 0.183 | 0.005 | 0.004 | 0.950 | 0.921 | 1.7 | 0.003 |
| 1k | R9.4.1 | even | 0.306 | 0.030 | 0.016 | 0.950 | 0.631 | 11.2 | 0.023 |
| 10k | R10.4.1 | even | 0.185 | 0.005 | 0.005 | 0.950 | 0.980 | 0.4 | 0.001 |
| 10k | R9.4.1 | even | 0.275 | 0.014 | 0.008 | 0.950 | 0.832 | 3.9 | 0.012 |
| 10k | R10.4.1 | log | 0.030 | 0.005 | 0.005 | 0.655 | 1.000 | 0.0 | 0.000 |
| 10k | R9.4.1 | log | 0.026 | 0.016 | 0.015 | 0.665 | 0.970 | 0.4 | 0.001 |
| 50k | R10.4.1 | even | 0.188 | 0.005 | 0.005 | 0.950 | 0.983 | 0.3 | 0.001 |
| 50k | R9.4.1 | even | 0.278 | 0.013 | 0.007 | 0.950 | 0.827 | 4.0 | 0.013 |

Species L1 on even communities is almost all in the two non-identifiable pairs (group
L1 0.005 on R10.4.1); recall 0.95 is one missed species in 20 (*S. argenteus*, Section
4.5). Log-staggered recall (0.61-0.67) is the detection limit (Section 5).

### 4.2 Database variant and -N (cell grid; `ablation_cell_db_N.tsv`, `fig_db_N.png`)

| Database / -N | Species L1 | Genus L1 | FP abundance |
|---|---|---|---|
| exact, 3 | 0.131 | 0.008 | 0.006 |
| exact, 10 | 0.143 | 0.012 | 0.002 |
| exact, 50 | 0.182 | 0.046 | 0.002 |
| c995, 3 | 0.520 | 0.317 | 0.156 |
| reps, 3 | 1.239 | 1.002 | 0.609 |

(Means over the four cell samples.) More secondary alignments do not help on the exact
database: species L1 and the number of false-positive species (49, 83, 174 summed over
the four samples) rise with -N. The -N 3 result cannot be derived from -N 50 output:
minimap2 breaks ties differently and even the primary alignment changes.

The reduced databases fail through contamination, not through lost resolution
(`contamination_db_variants_*.tsv`). On reps, 27-68% of primary alignments land in the
wrong genus, and 2,729 of 2,741 such reads (even R10.4.1) hit a sequence that GTDB has
in more than one genus; 37-82% of the reps estimate is false-positive abundance on such
sequences (c995: 2-21%). The exact database is protected by its redundancy: the true
species' many copies fill the -N slots and outscore the contaminant copy.

### 4.3 Error profile, floor, estimator, pruning (ofat; `ablation_ofat_paired.tsv`)

- *Error profile* ('all' alignments vs 'primary', D3): no metric changes at 1k or 10k;
  at 50k species L1 +0.0008 (R10.4.1, 6/6 samples worse). Primary stays.
- *Output floor* (0, 1e-4, 1e-3): no change at 1k or 10k reads, because the pruning
  threshold max(10/R, floor) is already 1e-3 at 10k (Emu's rule: 1/R up to 1,000
  reads, 10/R above). At 50k the floor is binding: 1e-4 or 0 raises recall by 0.05-0.07
  on log-staggered samples, with 2-8 more false-positive species.
- *Estimator* (EM vs VB, D2): the same to three decimals in every metric.
- *Pruning off*: recall 1.0 with hundreds of false-positive species (153-700 at 50k).

### 4.4 Read score (M4)

The per-base score (matches and errors scored per aligned base) raises recall slightly
at 1k and 10k (10k even: 0.990 vs 0.950) at a small false-positive cost, and lowers
species L1 at 50k in all 12 samples (R10.4.1 −0.012, R9.4.1 −0.002). On the real mocks
it changes L1 by −0.004 to +0.012 and adds 0-2 false-positive species
(`real_arms_paired.tsv`).

### 4.5 Near-neighbour pairs (`pairs_ties.tsv`)

For each pair, reads were classified by whether their best scores tie between the two
members, favour the true member or favour the other.

| Pair | Tied reads (R10.4.1 / R9.4.1, 10k even) | VB share of member A (realised 0.50) |
|---|---|---|
| *S. aureus* / *S. argenteus* | 22.8% / 4.1% | 1.000 / 1.000 |
| *E. coli* / *E. fergusonii* | 12.7% / 0.8% | 0.766 / 0.839 |
| *B. subtilis* / *B. spizizenii* | 0.2% / 0.5% | 0.529 / 0.765 |
| *S. pneumoniae* / *S. mitis* | 0.4% / 0.9% | 0.590 / 0.691 |
| *K. pneumoniae* / *K. variicola* | 0.4% / 0.1% | 0.485 / 0.546 |

*S. argenteus* and *E. fergusonii* reads carry no evidence against their partner (16S
copies identical to the partner's), while the partner's own reads favour the partner;
VB and EM (identical here) then give the shared reads to the partner, and *S. argenteus*
is estimated at 0. With more secondary alignments the *E. coli* share rises further
(0.97 at -N 10, 0.99 at -N 50). The per-base score and pruning off do not change the
split. Pairs whose members have reads of their own split close to truth on R10.4.1;
on R9.4.1 the higher error rate moves 11-19% of those reads to the wrong member. More reads
do not help: at 50k reads (seeds 1-3) the error of VB's split, |share of member A −
realised share|, averaged over shapes and -N, is within 0.02 of its 10k value for
every pair and chemistry (*S. argenteus* pair 0.25 at 1k, 10k and 50k). The
e2e fixture's 51 tied *E. coli* / *S. flexneri* reads are the same case as
*E. fergusonii*.

## 5. Detection limit (M7, W6; `lod_detection.tsv`, `lod_summary.tsv`, `fig_lod.png`)

Detection by realised reads of a species, ofat samples, species outside the pairs:
detected in at least 90% of cases (and in every higher bin) from 1-2 reads at 1k reads
(the 1/R threshold), 11-20 reads at 10k (median share 0.0015-0.0017) and 51-100 reads at
50k (0.0019, the 1e-3 floor). Pair members need 21-50 reads at 10k. With pruning off,
1-2 reads suffice at every depth, at the false-positive cost of Section 4.3. Detection
is therefore set by the thresholds (10/R and the floor), not by the Bayesian model.

## 6. Real mock communities (M9, W4, W13)

**Data.** ONT open data zymo_16s_2025.09 (R10.4.1, SQK-MAB114.24; basecall BAMs
only): ZymoBIOMICS Microbial Community Standard plus four ATCC gDNA spike-ins (MSPlus,
6 replicates) and Gut Microbiome Standard D6331 (GutMock, 4 replicates), sup (primary)
and hac basecalls; the Emu paper's Zymo mock SRR10391201 (MinION R9.4) and gut mock
SRR14326315; the preprint's Zymo mock SRR23636353. Each subsampled to 10,000 reads
(`seqkit shuffle`, seed 20261022).

**Truth** (`benchmarks/truth/`, README there): the Zymo manuals' 16S-only composition
(D6300/D6305, D6331), the spike-ins as presence only, and Emu's imputed ONT truth
(Supplementary Table 4) for its gut mock, a weaker reference.

**Arms** (`real_arms_summary.tsv`, `real_arms_paired.tsv`, per species
`real_species.tsv`): reference; -N 50; c995 and reps; per-base score; EM; pruning
off; unknown component δ 0.02; identity ≥ 97%; Q10 and Q15; fastplong 0.7.0 +
minimap2 2.31; NanoVI on Emu's database (with and without δ 0.02 and identity ≥ 97%).

Results follow the simulations (`real_arms_paired.tsv`). EM = VB. -N 50 removes 1-3
false-positive species per sample; it lowers species L1 on the R9.4 Zymo mocks (−0.05,
−0.08), changes it by less than 0.003 on the ONT 2025.09 data, raises it on the Emu gut
mock (+0.04) and loses one gut species in five samples. c995 and reps add 6-25
false-positive species per sample and raise L1 by 0.08-0.44. Pruning off reports 94-409
species. The unknown component (δ 0.02) removes 1-4 false-positive species per sample with
L1 within 0.01, but reports 14-35% of the reads as unknown in mocks whose species are
all in the database (Section 9); an identity floor of 97% removes every read of the R9.4 mocks.
Per sample group, reference setting:

| Group | Species precision | Recall | Species L1 | Genus L1 | FP abundance |
|---|---|---|---|---|---|
| MSPlus sup (6) | 0.71 | 1.00 | 0.18 | 0.20 | 0.02 |
| GutMock sup (4) | 0.77 | 0.91 | 0.61 | 0.57 | 0.03 |
| SRR10391201 (R9.4) | 0.47 | 1.00 | 0.43 | 0.26 | 0.19 |
| SRR23636353 | 0.62 | 1.00 | 0.35 | 0.31 | 0.12 |
| SRR14326315 | 0.59 | 1.00 | 0.28 | 0.16 | 0.12 |

Three findings explain most of the error:

- *The GutMock reads deviate from the theoretical composition.* Genus L1 is 0.57 with
  genus precision and recall 1.0; counting each read's best hit gives the same genus
  shares as VB (*Veillonella* 0.24, *Bacteroides* 0.21 against 0.16 and 0.10 expected;
  *Bifidobacterium* 0.02, *Limosilactobacillus* 0.01 against 0.09 and 0.10), the
  pattern of amplification bias. The comparators on the same reads (Section 10) show
  whether any method recovers the theoretical shares.
- *The Zymo standards' Bacillus is B. spizizenii in GTDB.* Zymo's reference genome of
  its "*B. subtilis*" strain is 99.97% identical over its whole length to the
  *B. spizizenii* type strain ATCC 6633 (GTDB's representative) and 87.8% identical over
  94% to *B. subtilis* 168; five of its ten 16S copies are identical only to
  *B. spizizenii* database sequences, the others to sequences shared with *B. subtilis*
  and relatives. The truth tables use *B. spizizenii* (`benchmarks/truth/README.md`).
  On R10.4.1 sup reads NanoVI reports it as *B. spizizenii* (MSPlus mean 0.144 of 0.174
  expected), with *B. subtilis* 0.016 and three other *Bacillus* species at about 0.003
  each as false positives. On the R9.4 mocks the reads spread over seven species of the
  complex (SRR10391201: *B. subtilis* 0.110, *B. spizizenii* 0.040, *B. halotolerans*
  0.032; SRR23636353, the preprint's mock: *B. spizizenii* 0.115, *B. halotolerans*
  0.053, *B. subtilis* 0.044), the error-rate effect seen in the simulated pair
  (Section 4.5). (W13.)
- *A contaminant sequence absorbed correctly.* In GutMock rep1, 297 reads have their
  primary alignment on *Staphylococcus epidermidis*; 278 of them on one GTDB sequence
  (43435) that is identical to a *Veillonella* 16S (a cross-genus identical sequence),
  and 305 of their secondary alignments score the same on *Veillonella*. Counting
  primary hits would report 3% *S. epidermidis*; VB reports none.

## 7. Read filtering and tool versions (M8, D5, D16)

fastplong's `--mean_qual` averages the Phred values arithmetically. Badread's qualities
are high (R9.4.1-like median Q23, 5th percentile 19.5; R10.4.1-like Q33 and 31.5), so Q10
removes nothing from the simulations; the simulated arms use Q20/Q22 (R9.4.1) and
Q32/Q33 (R10.4.1) instead, and the real data Q10 and Q15. Quality-filtered alignments
are taken from the reference alignment of the reads that pass (`harness.subset_bam`;
minimap2 maps each read independently, checked identical to realigning).

- *Simulations* (`quality_sims.tsv`, 40 samples at 10k): on R9.4.1 even communities,
  Q22 keeps 71% of the reads and raises species precision from 0.83 to 0.94
  (FP abundance 0.012 to 0.007); on R10.4.1 Q33 keeps 42%. Genus L1 rises (R10.4.1 even
  0.005 to 0.035 at Q33): the filter removes reads unevenly across species.
- *Real mocks*: Q10 and Q15 change species L1 by less than 0.008 and the false-positive
  count by at most one species; Q15 removes 1-16% of the reads.
- *fastplong 0.7.0 + minimap2 2.31* (D16): on every real sample, the same species and
  L1 within 0.0004 of the reference; assigned reads differ by at most 3 per 10,000.

## 8. Credible intervals (M2, D7, W2)

`intervals.py` computes three 95% intervals on the reference setting's scores for all
92 ofat samples: the Phase 5 mean-field interval (Beta from the VB posterior), a read
bootstrap (100 resamples, VB rerun on each) and a collapsed Gibbs sampler on the same
sparse likelihood (an exact-posterior reference). Coverage is the share of
(sample, species) pairs whose interval contains the design abundance, renormalised over
the reported species (`intervals_coverage*.tsv`, `fig_intervals.png`).

| Method | 1k | 10k | 50k | 10k, pair members |
|---|---|---|---|---|
| mean-field | 0.93 | 0.81 | 0.38 | 0.42 |
| bootstrap | 0.91 | 0.79 | 0.39 | 0.41 |
| collapsed Gibbs | 0.93 | 0.81 | 0.39 | 0.42 |

(Species outside the pairs, except the last column.) The three methods agree in width
and coverage at every depth, so the mean-field interval is not too narrow relative to
the exact posterior. Coverage falls with depth because the intervals shrink as 1/√R
while the estimates keep a systematic error from read assignment (the median error is
0.4 half-widths at 1k, 0.6 at 10k and 1.2 at 50k for species outside the pairs, 3.3 at
50k for pair members). The intervals describe sampling uncertainty only.

`shared_read_fraction` flags the misses: AUC 0.75 (1k), 0.84 (10k), 0.79 (50k)
(`intervals_warning_auc.tsv`); above 0.1, coverage is at most 0.6 at 1k and 0.13 at 10k
(`intervals_coverage_by_shared_read_fraction.tsv`). Median runtime per sample at 1k / 10k /
50k reads: mean-field 0.02 / 0.09 / 0.7 s, Gibbs 0.04 / 0.3 / 1.2 s, bootstrap 1.6 / 8.9 /
80 s (`intervals_runtime.tsv`; loaded host).

## 9. Unknown-taxon component and filters (M5, D10)

Design: `docs/design/unknown_component.md`. Clade exclusion (`clade_exclusion.py`):
each of the 20 community species removed in turn from the exact database (all its
sequences), the six 10k even samples (seeds 1-3, both chemistries) realigned to each
reduced database (120 alignments), then the excluded species' reads followed through
14 arms: none; identity ≥ 85/90/95/97/98/99%; aligned length ≥ 1,000/1,200/1,400; the
unknown component with δ 0.01/0.02/0.03/0.05 (`clade_exclusion_summary.tsv`,
`fig_clade.png`).

Without filtering, the excluded species' reads go to species absent from the sample
(58-66% of them) or, when the pair partner is present, mostly to the partner (49-57%). None of the three
approaches separates them from reads of species in the database: wherever an arm
removes at least 1% of the reads, it removes the excluded species' reads at 0.9-1.8
times the rate at which it removes reads of in-database species.

| Arm | Excluded reads removed (R10.4.1 / R9.4.1) | In-database reads removed |
|---|---|---|
| unknown, δ 0.01 | 0.22 / 0.42 | 0.17 / 0.39 |
| unknown, δ 0.02 | 0.08 / 0.31 | 0.06 / 0.28 |
| identity ≥ 99% | 0.64 / 1.00 | 0.55 / 1.00 |
| identity ≥ 95% | 0.02 / 0.62 | 0.02 / 0.60 |
| length ≥ 1,400 | 0.02 / 0.05 | 0.02 / 0.05 |

Identity filters are chemistry-dependent: 97% removes 93% of all R9.4.1 reads (and
every read of the R9.4 published mocks); only 99% does anything on R10.4.1. The
design is the hardest case: every excluded species has GTDB relatives within the
read error of its 16S, so per-read evidence cannot tell a missing species from a
present relative. A genus-level exclusion would be an easier test; it was not run.

On the real mocks, where every species is in the database, δ 0.02 reports 14-20% of
the reads as unknown on the ONT 2025.09 sup data, 23% on SRR23636353 and 34-35% on the
R9.4 mocks (`real_nanovi.tsv`, arm unk0.02), against 6% (R10.4.1) and 28% (R9.4.1) of
in-database reads on the simulations: real reads have heavier error tails than
Badread's, and the floor, set from each sample's mean error, sits inside them. The
database species, renormalised over themselves, keep their L1 (within 0.01), and 1-4
false-positive species per sample disappear. On Emu's database the same arm loses
*Dorea longicatena* on the Emu gut mock and raises species L1 by 0.09-0.26. Reported
as a share of the sample, the unknown row would read as 14-35% novel taxa in a mock
with none.

## 10. Comparators (D8, W4, W5)

Source data: `comparators_per_sample.tsv`, `comparators_summary.tsv`,
`comparators_paired.tsv`, `fig_comparators.png`. Every tool gets the same reads
(Emu and Kraken2 the reads NanoVI aligns; wf-16s the unfiltered reads, which it filters
itself), 10,000 per sample. Names from NCBI databases are mapped to GTDB r226 before
scoring, through the species' genomes in GTDB's metadata and, for older names such as
Emu's *Ruminococcus gnavus* and *Paeniclostridium sordellii*, NCBI synonyms
(`benchmarks/names.py`). Species-level scores across databases are not fair: one strain can carry
different names (Zymo's *Bacillus* is *B. subtilis* in NCBI and *B. spizizenii* in
GTDB, 17% of the MSPlus mock). The comparison therefore uses the identical-16S group
level (species whose 16S sequences are identical in GTDB merged) and reports species
level alongside. Every tool except wf-16s in minimap2 mode already drops species
below 0.1% at this depth; wf-16s minimap2 reports a tail (25.5 false-positive species
per real sample, 13.6 at 0.1% or more; `fp_count_ge_0.1pct`).

Real mocks (13 samples; Emu on NanoVI's database on 5):

| Tool / database | Species L1 | Group L1 | Group recall | False-positive groups | False-positive abundance |
|---|---|---|---|---|---|
| NanoVI (VB) / exact GTDB | 0.352 | 0.333 | 0.972 | 3.6 | 0.028 |
| NanoVI (EM) / exact GTDB | 0.353 | 0.333 | 0.972 | 3.8 | 0.028 |
| NanoVI (VB) / Emu's database | 0.569 | 0.413 | 0.967 | 2.9 | 0.104 |
| Emu 3.6.2 / Emu's database | 0.568 | 0.415 | 0.967 | 1.2 | 0.096 |
| Emu 3.6.2 / exact GTDB (5) | 0.355 | 0.321 | 0.964 | 2.6 | 0.033 |
| Kraken2 + Bracken / ncbi_16s_18s | 0.824 | 0.811 | 0.862 | 25.0 | 0.349 |
| wf-16s v1.6.1 minimap2 / ncbi_16s_18s | 0.711 | 0.494 | 0.934 | 16.0 | 0.200 |
| wf-16s v1.6.1 kraken2 / ncbi_16s_18s | 0.831 | 0.818 | 0.852 | 24.9 | 0.355 |

Simulations (12 samples, 10k reads, seeds 1-3; Emu on the 4 seed-1 samples):

| Tool / database | Species L1 | Group L1 | Group recall | False-positive groups |
|---|---|---|---|---|
| NanoVI (VB) / exact GTDB | 0.131 | 0.010 | 0.869 | 0.5 |
| NanoVI (EM) / exact GTDB | 0.131 | 0.010 | 0.869 | 0.4 |
| Emu 3.6.2 / Emu's database (4) | 0.222 | 0.011 | 0.845 | 0 |
| Emu 3.6.2 / exact GTDB (4) | 0.229 | 0.088 | 0.845 | 0 |
| NanoVI (VB) / Emu's database (4) | 0.202 | 0.010 | 0.845 | 1.0 |
| Kraken2 + Bracken / ncbi_16s_18s | 1.130 | 1.049 | 0.805 | 72.1 |
| wf-16s v1.6.1 minimap2 / ncbi_16s_18s | 0.846 | 0.177 | 1.000 | 6.2 |
| wf-16s v1.6.1 kraken2 / ncbi_16s_18s | 1.137 | 1.065 | 0.805 | 75.3 |

<!-- pending: NanoCLUST -->

Paired on the same samples (`comparators_paired.tsv`):

- **NanoVI as shipped vs Emu as shipped** (each on its own database). Real: group L1
  lower with NanoVI on all 13 samples (mean 0.083, range 0.020-0.201); Emu reports 2.5
  fewer false-positive groups (fewer on 12 of 13); group recall equal on 12 of 13 (on
  the Emu-paper gut mock Emu misses *Romboutsia hominis*, which its database lacks; it
  reports *R. timonensis*). Simulations: group L1 equal (±0.016), species L1
  lower with NanoVI on all 4 (mean 0.089), species recall 0.10 higher with NanoVI.
- **Same database, both ways.** On Emu's database NanoVI and Emu give the same group L1
  (within 0.013 on every sample) and the same recall; NanoVI reports 1.8 more
  false-positive groups (more on all 13). On the simulations, on Emu's database,
  group L1 is equal (within 0.018), species L1 is lower with NanoVI on all 4 (mean
  0.021) and NanoVI reports 1 more false-positive group. On NanoVI's database (5 real samples) group
  L1 differs by -0.011 to +0.037 (mean 0.009 in NanoVI's favour) and Emu reports 1.4
  fewer false-positive groups; on the simulations Emu's group L1 is 0.076 higher. The
  difference on the simulations comes from Emu's -N 50: with -N 50 the GTDB
  contaminant copies (Section 4.2) get alignment slots and pull *Enterococcus* and
  *Acinetobacter* reads into the *Klebsiella* and *E. coli* groups, as NanoVI's own
  -N 50 cell does (group L1 0.046 vs 0.012 at -N 3).
- **VB vs EM in NanoVI's code**: no difference (group L1 within 0.002, 0-1 false-positive
  group).
- **Kraken2 + Bracken**: worse than NanoVI on every sample (group L1 +0.48 real, +1.04
  simulated; 21 and 72 more false-positive groups). wf-16s in kraken2 mode gives the
  same results (it runs the same tools on the same database).
- **wf-16s, minimap2 mode** (best hit per read, identity ≥ 95%, coverage ≥ 90%): group
  L1 higher than NanoVI's on 10 of 13 real mocks (mean +0.16) and on all 12
  simulations (+0.17), with 12 and 6 more false-positive groups; it detects more of
  the rare simulated species (group recall 1.00 vs 0.87), since it applies no
  abundance threshold.

What this supports: NanoVI's lower error than Emu-as-shipped comes from the database
(exact GTDB r226 vs Emu's rrnDB + NCBI RefSeq), not from the estimator. With the same
database the two give the same abundances, and Emu reports fewer false positives.
That difference comes from the alignment settings (Emu: -N 50, k 15; NanoVI: -N 3,
k 21), not from the estimator: given the same SAM, Emu, NanoVI's VB and NanoVI's EM
report the same species and the same false positives (`timing_matched_accuracy.tsv`:
2 samples × 4 settings of -N and k on the exact database; identical false-positive
species and groups in all 8, species L1 within 0.0009, group L1 within 0.0002).

## 11. Runtime (W3)

Source data: `timing_steps.tsv` (per step, medians of 3 repeats with min and max),
`timing_tools.tsv` (per tool), `fig_timing.png`; script `benchmarks/timing.py`. Two
real samples of about 10,000 filtered reads (ONT 2025.09 MSPlus sup rep 1; Emu
paper SRR10391201, R9.4.1), each aligned once per setting (-N 3 or 50, k 15 or 21) to
the exact database, and the same SAM given to NanoVI (BAM conversion, likelihood, VB
or EM) and to Emu 3.6.2 (`emu abundance` on the SAM). Emu as shipped is timed
separately: its minimap2 2.30 from the FASTA (it builds the index every run), -N 50,
k 15, its own database. NanoVI as shipped is the -N 3, k 21 row with a cached index;
the index build (25 CPU s, 1.7 GB, once per k) is reported apart. Every step ran on
4 pinned CPUs. The host was shared (1-minute load 35-53 on 80 threads), so CPU time
(user + system) is the number to compare; wall times are in the tables with the load.

| Sample | NanoVI as shipped (-N 3, k 21) | Emu as shipped | Emu's settings (-N 50, k 15), exact database: NanoVI / Emu |
|---|---|---|---|
| SRR10391201 (R9.4.1) | 870 CPU s (246 s wall) | 2,750 CPU s (816 s wall) | 3,120 / 3,160 CPU s |
| MSPlus sup rep 1 (R10.4.1) | 5,097 CPU s (1,425 s wall) | 7,077 CPU s (1,845 s wall) | 9,822 / 9,849 CPU s |

- **The alignment is the cost**: 96-99.9% of CPU time for both tools at every
  setting. k 15 costs 1.9-3.2 times the CPU of k 21; -N 50 adds 1-60% over -N 3.
  Alignment cost depends on the sample (0.09 vs 0.52 CPU s per read at -N 3, k 21).
  Peak memory is the alignment's (0.9-1.6 GB with the exact database's index).
- **After the alignment** NanoVI (BAM conversion, likelihood, VB) takes 4.5-5.5 CPU s at
  -N 3 and 13-19 at -N 50; Emu on the same SAM 8.6-11.8 and 41-59 (1.9-3.1 times
  NanoVI). NanoVI's EM adds about 2 CPU s over VB. Peak memory below 120 MB for NanoVI
  and 190 MB for Emu.
- **As shipped**, NanoVI used 1.4-3.2 times less CPU than Emu on these samples, with
  a database 4.7 times larger (232,447 vs 49,243 sequences). The saving comes
  from -N 3 and k 21, not from the estimator; at Emu's settings on the exact database
  the two cost the same.
- Accuracy at these settings: §10 (identical between the estimators on the same SAM).
  Changing -N and k moves the false-positive count by 1-2 groups on these two samples,
  in both directions (`timing_matched_accuracy.tsv`); -N's effect on the simulations is
  in §4.

## 12. Limitations of this benchmark

- One simulated community; its two non-identifiable pairs dominate species-level error
  on even communities.
- Badread's quality scores are higher than real basecalls', so the quality-filter
  thresholds on simulations are not comparable to Q10/Q15 on real reads.
- Real data subsampled to 10,000 reads; the cell grid on four samples; 50k reads with
  three seeds.
- Wall times from the grid were measured on a loaded shared host and are not used;
  Section 11 uses CPU time from pinned runs.
- The clade exclusion removes species, not genera.
- Simulated depths and shapes with the same seed share each species' reads (Section 3),
  so differences between depths are paired, not independent, comparisons.
