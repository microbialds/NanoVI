# NanoVI on clinical 16S cohorts (Phase 9, October 2026)

This document records how the clinical rerun of the manuscript (decision D34) was run
and what it found. With no ground truth in clinical samples, it measures concordance:
between NanoVI and Emu on the same reads, with each database and with the same
alignments, and with the profiles the cohorts' authors published. Every table cited
here is in `docs/benchmarks/data/` (`clinical_*`), the figure is the manuscript's
Figure 8 (`docs/paper/make_figures.py`, `fig8_clinical`), and the code is
`benchmarks/clinical.py`. Data, runs and
per-sample outputs are under `~/pipelines/NanoVI_benchmarks/phase9` (D11) and are not
in the repository.

## 1. Cohorts

| Key | BioProject | Samples | Reads | Published reference | Labels |
|---|---|---|---|---|---|
| vaginal_cl | PRJNA1232418 (Oliva-Arancibia et al. 2025, *Sci Rep*; the preprint's cohort) | 20 vaginal self-samples (M1-M20), MinION | 5,000-42,993 per sample, all used | Genus and species tables of the paper (Supplementary Tables 4 and 5, Emu on Emu's database) | none |
| vaginal_bv | PRJNA723982 (Curry et al. 2022, *Nat Methods*; Emu's clinical example) | 12 vaginal swabs, 6 controls (P01-P06) and 6 with bacterial vaginosis (P07-P12), MinION 16S | 10,000 per sample (subsampled) | Emu's per-sample tables and the Kraken2 + Bracken profiles of the same swabs' ONT shotgun reads (gitlab.com/treangenlab/emu-benchmark, `results_clinicial_application/`) | BV / control |
| crc | PRJNA911189 (Aja-Macaya et al. 2025, *Sci Rep*) | 117 stool samples, 89 colorectal cancer (CCR) and 28 healthy (SANO), MinION R10.4.1, Dorado sup (the ONT 16S runs in the BioProject; the paper reports 123 subjects) | 10,000 per sample (subsampled) | The four species of the paper's smallest CRC model (its Table 1: *F. nucleatum*, *P. micra*, *B. fragilis*, *A. butyriciproducens*) | CRC / healthy |

Only the ONT 16S reads are used (no Illumina, decision D34). Subsampling: `seqkit shuffle
-s 20261022` then the first 10,000 reads, as for the real mocks (`methods_2026.md` §6).
Manifests (ENA read-run reports with MD5s) are kept beside the data.

## 2. Arms

| Arm | What runs | Database |
|---|---|---|
| `nanovi_gtdb` | `nextflow run main.nf -profile docker`, every default (NanoVI as shipped) | GTDB r226, default build (suffixes kept; `bin/build_gtdb_db.py`) |
| `nanovi_emudb` | the same pipeline, `--db` Emu's database | Emu's default database (rrnDB v5.6 + NCBI 16S RefSeq) |
| `emu` | Emu 3.6.2 as shipped (`emu abundance`, -N 50, k 15, its minimap2 2.30), on the reads NanoVI kept after length filtering | Emu's default database |
| `emu_same_sam` | Emu 3.6.2 on the alignments of `nanovi_emudb` (BAM converted to SAM), so both tools score the same alignments | Emu's default database |

`nanovi_emudb` against `emu_same_sam` isolates the estimator; `nanovi_emudb` against `emu`
adds the alignment settings; `nanovi_gtdb` against `emu` is the two tools as shipped.

## 3. Metrics

- **Bray-Curtis dissimilarity** between two profiles of the same sample (0 identical, 1
  disjoint), at species and genus level, over the assigned reads (Unassigned and Below
  detection rows left out, each profile renormalised).
- **Names across databases.** NCBI names (Emu's database, the published tables, Bracken)
  are mapped to GTDB r226 species with the Phase 7a crosswalk (`benchmarks/names.py`:
  the GTDB species most of the NCBI species' genomes are placed in, plus NCBI synonyms);
  GTDB names are compared with suffixes stripped. Genus comparisons across databases use
  the genus of the mapped GTDB species (the first word of its name); within Emu's database
  the genus is Emu's taxonomy's, or the first word of the name for species it gives no
  genus. The abundance that cannot be mapped is reported per sample (`unmapped`) and left
  out of the NCBI-named profile, which is renormalised; the GTDB profile it is compared
  with keeps those reads under the GTDB species they went to, so a large unmapped share
  adds to the dissimilarity.
- **Dominant taxon** per sample (the most abundant genus) and its agreement between arms.
- **Labels.** vaginal_bv: the *Lactobacillus* share of each profile in control and BV
  swabs, its range per group and the AUC of a Mann-Whitney U test (control above BV),
  for every 16S arm, Emu's published tables and the shotgun Bracken genera
  (`clinical_bv_lactobacillus.tsv`; bacterial genera only, by their lineage in the NCBI
  taxdump of 2025-01-01, as the shotgun reads also hit human and other eukaryotic genera,
  up to 14% of a swab's non-human genus-level reads). crc: for each arm, a two-sided
  Mann-Whitney U test per species present in at least 10% of samples with
  Benjamini-Hochberg correction (`clinical_crc_species.tsv`); the counts below q 0.05 and
  0.10 per arm (`clinical_crc_tests.tsv`); and the direction and q of the paper's four
  species (`clinical_crc_paper_species.tsv`, NCBI names mapped as above). GTDB r226
  splits two of them, so each is tested as the GTDB species holding at least 5% of its
  NCBI genomes, merged before the tests: *F. nucleatum* (with *F. polymorphum*, *F.
  animalis* and *F. vincentii*) and *B. fragilis* (with *B. hominis*).
- **Database content.** Each taxon's mean contribution to the Bray-Curtis dissimilarity
  between NanoVI on GTDB and NanoVI on Emu's database (GTDB names, the ten largest per
  cohort and level), with the number of species of its genus in each database
  (`clinical_db_contribution.tsv`).

## 4. Results

Medians over samples from `clinical_summary.tsv` (BC = Bray-Curtis dissimilarity; max in
brackets); figure: `docs/paper` Figure 8 (`fig8_clinical`).

| Comparison | vaginal_cl (20) | vaginal_bv (12) | crc (117) |
|---|---|---|---|
| NanoVI vs Emu, same alignments: species BC | 0.0040 (0.020) | 0.0001 (0.001) | 0.0014 (0.006) |
| same: dominant genus agrees | 100% | 100% | 100% |
| NanoVI vs Emu, each as shipped, Emu's database: species BC | 0.029 (0.18) | 0.0054 (0.12) | 0.067 (0.40) |
| same: genus BC | 0.0035 | 0.0003 | 0.044 |
| same: dominant genus agrees | 95% | 100% | 93% |
| Emu 3.6.2 rerun vs the authors' Emu tables: species BC | 0.0061 (0.17) | 0.0067 (0.036) | — |
| NanoVI (Emu's database) vs the authors' Emu tables: species BC | 0.029 (0.21) | 0.0098 (0.12) | — |
| NanoVI, GTDB r226 vs Emu's database: species BC | 0.078 (0.80) | 0.070 (0.77) | 0.48 (0.80) |
| same: genus BC | 0.016 | 0.0041 | 0.24 |
| same: dominant genus agrees | 85% | 83% | 81% |
| NanoVI on GTDB r226 vs Emu as shipped: species BC | 0.079 | 0.070 | 0.48 |

**Estimator.** On the same alignments NanoVI and Emu gave the same profiles (species BC at
most 0.020 in any sample, the same dominant genus in all 149), as on the mock communities
(`methods_2026.md` §10).

**Alignment settings.** Each tool as shipped on Emu's database (NanoVI -N 3 -k 21, Emu
-N 50 -k 15) differed more, most in stool (median species BC 0.067). Of the 8 stool
samples with a different dominant genus, 4 are *Butyricimonas* (NanoVI) against reads Emu
gave to unnamed Bacteroidales species (no genus in Emu's taxonomy), and
one (SRR30050069) has 17% *Bacillus thuringiensis* in NanoVI's profile where Emu as
shipped reports *Faecalibacterium* first (`clinical_pairs.tsv`, its rel-abundance tables).
Our Emu rerun reproduces the authors' published Emu tables (median species BC 0.006 in
both vaginal cohorts; the published tables give their main taxa plus 'Others').

**Database.** NanoVI on GTDB r226 differed from both Emu-database arms by the same amount
(species BC 0.48 against NanoVI on Emu's database and against Emu, in stool), so the
difference is the database, not the tool. It comes from what the databases contain
(`clinical_db_contribution.tsv`): species GTDB has and Emu's database lacks
(*Amygdalobacter indicium*, BVAB2, 26% of the genus-level difference in vaginal_cl; 0
species of the genus in Emu's database), genera GTDB splits into many species
(*Faecalibacterium*: 32 species in GTDB, 1 in Emu's database; the top species-level
contributor in stool), placeholder species and genera (*Dialister* sp001553355,
*Cryptobacteroides*, CAG-273) and species GTDB separates (*L. gasseri* and *L.
paragasseri*). The genus-level difference is smaller in the vaginal cohorts (median 0.004
-0.016) than in stool (0.24), where it is spread over many genera (the ten largest
contribute 30%). Part of Emu's database profile has no GTDB name in the crosswalk and is
left out of that profile (median 0-0.5% in the vaginal cohorts and 7.5-8.7% in stool, up
to 56% in P11, *Peptostreptococcaceae bacterium oral taxon 929*; `clinical_pairs.tsv`,
`unmapped_b`); NanoVI's GTDB profile keeps those reads, so they add to the dissimilarity.
The genus counts compare each database's own genus names, so a genus GTDB renamed or
split (*Lactobacillus*: 68 species in GTDB, 192 in Emu's database) is counted under different concepts.

**Bacterial vaginosis.** Every 16S arm separated the 6 control from the 6 BV swabs by
*Lactobacillus* share (controls 0.986-1.000, BV 0.000-0.967; AUC 1.0, Mann-Whitney p
0.005), as did Emu's published tables; the shotgun Bracken genera did not quite (controls
0.602-0.982, BV 0.019-0.645; AUC 0.97) (`clinical_bv_lactobacillus.tsv`). The shotgun
reference is shallow (170-14,290 bacterial reads per swab, `clinical_samples.tsv`) and
names *Gardnerella* first in three BV swabs where the 16S profiles name *Megasphaera* or
*Lactobacillus*; NanoVI and Emu are equally far from it (genus BC vs shotgun within 0.002
of each other in every swab).

**Colorectal cancer.** No species passed q < 0.05 with NanoVI on either database or with
Emu on NanoVI's alignments; Emu as shipped had 4 (q 0.041), with 167-198 species tested
per arm (`clinical_crc_tests.tsv`). *Gemmiger variabilis* (lower in CRC) was the lowest-q
species on Emu's database with every arm (q 0.041-0.069). Of the paper's four species,
*P. micra* (AUC 0.60-0.61, absent from healthy samples) and *B. fragilis* (0.63-0.64) were
higher and *A. butyriciproducens* (0.33-0.38) lower in CRC with every arm and both
databases, as in the paper, though not significant after correction (q 0.11-0.41);
*F. nucleatum* was in at most 2% of samples and was not tested
(`clinical_crc_paper_species.tsv`). With 10,000 reads per sample and 28 healthy samples
this is a test of concordance between the arms, not a replication of the paper.
