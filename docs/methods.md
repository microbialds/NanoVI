# NanoVI methods

How NanoVI turns full-length 16S Nanopore reads into species abundances, what it
takes from Emu and what it changes, the evidence behind each default, how the GTDB
reference is built, and what the method cannot do. Numbers come from the 2026
benchmark, [benchmarks/methods_2026.md](benchmarks/methods_2026.md) (cited as §), and
its source data in `benchmarks/data/`. Decisions are cited by their roadmap keys
(D2-D31). The output files are described in [output.md](output.md).

**In one paragraph.** NanoVI scores each read against the species it aligns to with
Emu's CIGAR-based read score, and estimates species abundances with a Bayesian
(Dirichlet-categorical) model fitted by mean-field variational inference, with Emu's
pruning of low-abundance species. Its point estimates match those of Emu's EM: the two
estimators agreed to three decimals on the simulations and within 0.002 in group L1 on
the real samples. With its exact
GTDB r226 reference, NanoVI was more accurate than Emu with Emu's default database on
all 13 real mock communities tested; with the same reference the two agree, and the
gain comes from the database, not the model (§10). The posterior gives sampling
intervals and a per-species ambiguity signal (`shared_read_fraction`), not calibrated
uncertainty (§8). Species whose 16S genes are identical cannot be separated; the
output names them.

## 1. Pipeline

| Step | Tool | What it does (defaults) |
|---|---|---|
| FILTER_READS | fastplong 0.7.0 | Keeps reads of 500-2,000 bp (`--min_length`, `--max_length`). The quality filter is off (`--disable_quality_filtering`); `--min_quality Q` turns on fastplong's mean-quality filter (`--mean_qual Q`) |
| INDEX_DB | minimap2 2.31 | Indexes the database's `species_taxid.fasta` with `-x map-ont -k 21`, once; the index is cached under a name that records the FASTA's MD5, the preset, k and the minimap2 version |
| ALIGN_READS | minimap2 2.31, samtools 1.24 | `minimap2 -ax map-ont -N 3 -p 0.9 -K 500M -k 21`: the primary alignment and up to 3 secondary alignments scoring at least 90% of it; BAM output |
| LIKELIHOOD | `bin/read_likelihood.py` | The sample's error profile and each read's score against each candidate species (section 2) |
| RUN_VI | `bin/estimate_abundance.py` | Abundances by variational inference with pruning (sections 3-4); `--estimator em` runs Emu's EM instead (section 5) |
| WRITE_OUTPUT, GENERATE_PHYLOSEQ_TABLES | `bin/write_abundance_table.py`, `bin/generate_phyloseq_tables.py` | Per-sample tables with lineage, counts, `identical_16s_species` and `shared_read_fraction` (section 6); consolidated phyloseq tables |

## 2. Read score (Emu's)

The read score is Emu's (Curry et al., *Nature Methods* 2022), unchanged; the code is
adapted from Emu's (`bin/likelihood.py`, credited in the README and
`THIRD_PARTY_LICENSES`).

**Error profile.** Over the sample's primary mapped alignments (not secondary or
supplementary ones, as in Emu; D3), count insertions $n_I$, deletions $n_D$,
soft-clipped bases $n_S$ and mismatches $n_X = \mathrm{NM} - n_I - n_D$ (NM is the edit
distance minimap2 reports). Each operation's share is

$$
p_o = \frac{n_o}{n_I + n_D + n_S + n_X}, \qquad o \in \{I, D, S, X\}.
$$

An operation with no count in the sample is dropped from the profile, and an alignment
that uses it gets no score.

**Score of an alignment.** For alignment $a$ of read $r$, with operation counts
$n_o(a)$ and aligned length $L(a) = n_M(a) + n_I(a) + n_D(a) + n_S(a)$ (M, I, D and S
CIGAR operations),

$$
\ell(a) = \frac{L_{\max}(r)}{L(a)} \sum_{o} n_o(a) \, \ln p_o ,
$$

where $L_{\max}(r)$ is the longest aligned length among all of the read's alignments.
The read's score for species $s$ is its best alignment to any reference sequence of
$s$: $\ell_{rs} = \max_{a \in A(r,s)} \ell(a)$. Species the read does not align to have
no score (zero likelihood).

This is a pseudo-likelihood: matches are not scored, and the scaling by
$L_{\max}/L(a)$ penalises shorter alignments in proportion. A per-base error model with
matches was tested and not adopted (section 8, D29).

$R$ below is the number of reads with at least one scored alignment. Reads with no
scored alignment are reported as `Unassigned`.

## 3. Model and variational inference

**Model.** Let $S$ be the number of species with at least one scored read in the sample
and $\pi = (\pi_1, \dots, \pi_S)$ their relative abundances among the reads. With
$\alpha_0$ = `--alpha_prior` (default 1):

$$
\pi \sim \mathrm{Dirichlet}\!\left(\tfrac{\alpha_0}{S}, \dots, \tfrac{\alpha_0}{S}\right), \qquad
z_r \mid \pi \sim \mathrm{Categorical}(\pi), \qquad
p(x_r \mid z_r = s) \propto e^{\ell_{rs}} ,
$$

where $z_r$ is the species read $r$ came from. The prior adds $\alpha_0 = 1$
pseudo-read in total, whatever $S$.

**Mean-field updates.** NanoVI approximates the posterior by
$q(\pi, z) = q(\pi) \prod_r q(z_r)$, with $q(\pi) = \mathrm{Dirichlet}(\alpha)$ and
$q(z_r = s) = \phi_{rs}$, and alternates the coordinate-ascent (CAVI) updates

$$
\phi_{rs} = \frac{\exp\!\big(\ell_{rs} + \psi(\alpha_s) - \psi(\textstyle\sum_j \alpha_j)\big)}
{\sum_{t} \exp\!\big(\ell_{rt} + \psi(\alpha_t) - \psi(\textstyle\sum_j \alpha_j)\big)}, \qquad
\alpha_s = \frac{\alpha_0}{S} + \sum_r \phi_{rs},
$$

where $\psi$ is the digamma function and the sums over $t$ run over the read's
candidate species. The evidence lower bound

$$
\mathcal{L} = \sum_{r,s} \phi_{rs}\big(\ell_{rs} + \mathbb{E}_q[\ln \pi_s] - \ln \phi_{rs}\big)
+ \ln B(\alpha) - \ln B(\alpha_0/S) + \sum_s \big(\tfrac{\alpha_0}{S} - \alpha_s\big)\mathbb{E}_q[\ln \pi_s],
$$

with $\mathbb{E}_q[\ln \pi_s] = \psi(\alpha_s) - \psi(\sum_j \alpha_j)$ and
$\ln B(\alpha) = \sum_s \ln\Gamma(\alpha_s) - \ln\Gamma(\sum_s \alpha_s)$, is computed
at every iteration as a check (a fall of more than $10^{-9}|\mathcal{L}|$ prints a
warning); it does not decide convergence.

**Start and stop.** The first fit starts from $\alpha_s = \alpha_0/S + R/S$ (reads spread
evenly). The loop stops when the largest change in posterior mean abundance between
iterations, $\max_s |\Delta\, \alpha_s / \sum_j \alpha_j|$, is below `--vi_tolerance`
($10^{-7}$), or after `--vi_max_iterations` (1,000) iterations. A relative ELBO rule
was used before and stopped large samples after a couple of iterations (F05).

The likelihood is stored as padded reads × candidates arrays (a read's candidates come
from its primary, up to $N$ secondary and any supplementary alignments, so they are few), so memory grows with the number of reads and not with reads ×
species: 58 MB instead of 24 GB at 200,000 reads over 2,000 species, with identical
results (F23).

## 4. Pruning and reported abundances

**Pruning.** After each fit, species whose posterior mean $\mathbb{E}_q[\pi_s] =
\alpha_s / \sum_j \alpha_j$ is below

$$
\tau = \begin{cases} \max(1/R,\ m) & R \le 1{,}000 \\ \max(10/R,\ m) & R > 1{,}000 \end{cases}
$$

are removed, with $m$ = `--min_abundance` (0.001), and the rest are refitted with a new
prior $\alpha_0/S'$ over the $S'$ survivors, starting from their expected read counts.
This repeats until no species is below $\tau$. The 1/R and 10/R terms are Emu's rule;
with the default $m$ the threshold is 1/R up to 1,000 reads, 10/R up to 10,000 reads and
0.1% above, so nothing below 0.1% is reported at any depth. After 10 prunings the loop
fits once more and stops, with a warning if species remain below $\tau$; if every
species falls below $\tau$, the most abundant one is kept.

**Reported abundances.** The reported abundance of species $s$ is its posterior mean,
renormalised over the reported species. Reads whose candidate species were all pruned
(every $\phi_{rs}$ on a reported species is zero) are counted in the `Below detection`
row instead of being spread over the reported species. Estimated counts are
abundance × ($R$ − reads below detection), rounded by largest remainder so that they
add up exactly.

**The thresholds set the detection limit, not the model** (D26, §5). In the
simulations, species outside the near-neighbour pairs were detected in at least 90% of
cases from 1-2 reads at 1,000 reads (the 1/R threshold), 11-20 reads at 10,000 reads
and 51-100 reads at 50,000 reads (the 0.1% floor). With pruning off, 1-2 reads suffice
at every depth, at a cost of hundreds of false-positive species.

## 5. Relation to Emu's EM

Emu estimates $\pi$ by expectation-maximisation:

$$
\phi_{rs} = \frac{\pi_s e^{\ell_{rs}}}{\sum_t \pi_t e^{\ell_{rt}}}, \qquad
\pi_s = \frac{1}{R}\sum_r \phi_{rs}.
$$

The variational update replaces $\ln \pi_s$ by $\psi(\alpha_s) - \psi(\sum_j \alpha_j)$.
Since $e^{\psi(x)} \approx x - \tfrac12$ for $x$ above a few, this behaves like EM on
expected counts shifted by $\alpha_0/S - \tfrac12$: the two differ only for species with
a handful of expected reads, and pruning removes those (the threshold is at least one
read). In the benchmark the two gave the same abundances to three decimals on every
simulated sample and within 0.002 in identical-16S group L1 on every real sample, with
0-1 false-positive group apart; on the same alignments, NanoVI's VB, NanoVI's EM and Emu
3.6.2 reported the same false-positive species (§4.3, §10; D2). The Bayesian
formulation's contribution is the posterior used in section 6, not better point
estimates.

`--estimator em` runs Emu's EM loop (`bin/em_estimator.py`, adapted from Emu): start
from $\pi_s = 1/S$; inner loops of at most 20 iterations until
$\sum_s |\Delta \pi_s| < 5 \times 10^{-4}$, repeated until the log-likelihood rises by
less than `--vi_elbo_threshold` ($10^{-4}$); then one pruning below 1/R (10/R above
1,000 reads) and a refit, and, if `--min_abundance` is higher, a second cut at it and a
refit. It is kept as the baseline for the comparison above (D2).

## 6. Posterior summaries: `shared_read_fraction` and intervals

**`shared_read_fraction`** (always written, D27). From the final fit, over the
reported species $\mathcal{R}$, let $\tilde\phi_{rs} = \phi_{rs} / \sum_{t \in \mathcal{R}}
\phi_{rt}$. A read is shared for species $s$ when another reported species takes at
least 1% of it ($0 < \tilde\phi_{rs} \le 0.99$). Then

$$
f_s = \frac{\sum_r \phi_{rs}\, \mathbf{1}[0 < \tilde\phi_{rs} \le 0.99]}{\sum_r \phi_{rs}} ,
$$

the share of the species' expected reads that it splits with another reported species
(`shared_read_fraction` = $f_s$).
With `--estimator em` the responsibilities are EM's final $\phi_{rs}$ over the species
that survived pruning, so the column means the same. It is computed per row of the
table: when a row sums several tax IDs (`--binomial_names`, or a repeated name), $\tilde\phi$
in the indicator is the read's share on the whole row, so reads split only among those
tax IDs are not shared (RUN_VI reads the rows from `taxonomy.tsv`).
In the benchmark it flagged the estimates that missed the truth (AUC 0.75, 0.84 and
0.79 at 1,000, 10,000 and 50,000 reads); above 0.1, interval coverage was at most 0.6 at
1,000 reads and 0.13 at 10,000 (§8). The mean-field summaries cost 0.09 s per sample at
10,000 reads.

**Sampling intervals** (`--credible_intervals`, VB only, off by default; D7, D27). The
mean-field marginal of a reported species (or of a row $g$ that sums several tax IDs) is

$$
\pi_g \sim \mathrm{Beta}\big(\alpha_g,\ A - \alpha_g\big), \qquad A = \sum_{s \in \mathcal{R}} \alpha_s ,
$$

and `ci_lower`, `ci_upper` are its 2.5% and 97.5% quantiles. These describe the
sampling uncertainty of the read counts only. They shrink as $1/\sqrt{R}$, while the
error from assigning reads between close relatives does not, so their coverage of the
true abundance falls with depth:

| | 1,000 reads | 10,000 reads | 50,000 reads | 10,000, pair members |
|---|---|---|---|---|
| mean-field (shipped) | 0.93 | 0.81 | 0.38 | 0.42 |
| read bootstrap (100 resamples) | 0.90 | 0.79 | 0.39 | 0.41 |
| collapsed Gibbs (exact posterior) | 0.93 | 0.81 | 0.39 | 0.42 |

(95% intervals, species outside the near-neighbour pairs except the last column; §8.)
The three methods agree, so the mean-field interval is not too narrow for the model;
the model itself does not include the read-assignment error. They are therefore
reported as sampling intervals, opt-in, and not as calibrated uncertainty.

## 7. What differs from Emu, and why

| Component | Emu 3.6.2 | NanoVI | Why / evidence |
|---|---|---|---|
| Read filtering | none | length 500-2,000 bp; quality filter off | A mean-quality filter changed little on real data and distorted the composition where it removed reads (D28, §7) |
| Alignment | minimap2 `-ax map-ont -N 50 -p .9 -K 500M`, k 15, index built every run | `-N 3 -p 0.9 -K 500M`, k 21, index cached | On the simulations -N 3 gave the lowest species L1 on the exact GTDB database, with more false-positive species than -N 50 (1.6 against 0.25 per sample, false-positive abundance 0.005 against 0.001); on real mocks -N 50 removed 1-3 false-positive species per sample, with L1 lower on the R9.4 Zymo mocks and higher on the Emu gut mock; -N 3 is the cheapest (D21, §4.2, §6). k 15 costs 1.9-3.2 times the alignment CPU of k 21 (§11) |
| Error profile | primary mapped alignments | the same (D3) | Counting all alignments changed nothing at 1k/10k and was slightly worse at 50k (§4.3) |
| Read score | CIGAR pseudo-likelihood | the same (D29) | |
| Estimator | EM (point estimate) | mean-field VB with a Dirichlet prior; EM available | Same point estimates (§4.3, §10); the posterior gives section 6's summaries (D2) |
| Convergence | log-likelihood change | largest change in posterior mean < 1e-7 | Runs to the fixed point at every sample size (F05) |
| Pruning | 1/R, 10/R above 1,000 reads | the same, with a 0.1% floor (`--min_abundance`) | Floor kept: lower floors at 50k reads add recall 0.05-0.07 with 2-8 more false-positive species (D26, §4.3) |
| Reads that lose all species | not reported | `Below detection` row | Every read is accounted for (F06) |
| Extra outputs | | `shared_read_fraction`, `identical_16s_species`, optional sampling intervals, phyloseq tables | D24, D27 |
| Database | rrnDB + NCBI 16S RefSeq, NCBI taxonomy | GTDB r226 `ssu_all`, all distinct sequences per species, GTDB taxonomy (section 9) | The database is what makes NanoVI more accurate than Emu as shipped (§10) |

**Accuracy against Emu and other tools** (real mocks, 10,000 reads, 13 samples;
identical-16S group L1 / false-positive groups per sample; §10): NanoVI on its GTDB
database 0.333 / 3.6; Emu on its database 0.415 / 1.2; Emu on NanoVI's database
(5 samples) 0.321 / 2.6; NanoVI on Emu's database 0.413 / 2.9; EPI2ME wf-16s in minimap2
mode 0.494 / 16.0; NanoCLUST (12 Zymo mocks) 0.543 / 2.9; Kraken2 + Bracken 0.811 / 25.0.
NanoVI as shipped had lower group L1 than Emu as shipped on all 13 samples (mean
difference 0.083). On the same database the two give the same L1, and Emu reports 1-2
fewer false-positive groups; given the same alignments they report the same false
positives, so that difference comes from Emu's alignment settings (-N 50, k 15).

**Runtime** (§11). Alignment is 96-99.9% of the CPU time for both tools. As shipped,
NanoVI used 1.4-3.2 times less CPU than Emu on two real samples (870 vs 2,750 and 5,097
vs 7,077 CPU s), with a database 4.7 times larger; at Emu's settings on the same
database they cost the same. The saving comes from -N 3 and k 21.

## 8. Defaults and their evidence

Every default below was kept or set after the 2026 benchmark (Phase 7a) and confirmed
on the released code (Phase 7b, §13).

| Setting | Default | Decision | Evidence |
|---|---|---|---|
| Quality filter | off (length only) | D5, D28 | Q10 and Q15 changed real-mock species L1 by less than 0.008 and the false-positive count by at most one species; on the simulations, thresholds high enough to matter removed reads unevenly across species (R10.4.1 even, genus L1 0.005 → 0.035 at Q33) (§7) |
| Tool versions | fastplong 0.7.0, minimap2 2.31 | D30 | Against fastplong 0.2.2 and minimap2 2.24: the same species on every real sample, species L1 within 0.0004, assigned reads within 3 per 10,000 (§7, §13) |
| Secondary alignments | `-N 3` | D21 | Simulations: species L1 0.130 / 0.143 / 0.182 at -N 3 / 10 / 50 on the exact database (genus L1 0.008 / 0.012 / 0.046), while false-positive species fall (1.6 / 0.5 / 0.25 per sample): at -N 50 GTDB's cross-genus identical copies take alignment slots and move reads between species of the community (§4.2). Real mocks: -N 50 removes 1-3 false-positive species per sample, lowers species L1 on the R9.4 Zymo mocks (−0.05, −0.08), changes it by less than 0.003 on the ONT 2025.09 data, raises it on the Emu gut mock (+0.04) and loses one gut species in five samples (§6). -N 3 is also the cheapest |
| Error profile | primary alignments | D3, D31 | 'All alignments' changed no metric at 1k or 10k reads; species L1 +0.0008 at 50k (§4.3) |
| Read score | Emu's | D29 | A per-base score raised recall on the simulations (10k even: 0.990 vs 0.950) but added false positives on real mocks (L1 −0.004 to +0.012, 0-2 more false-positive species) (§4.4) |
| Identity and aligned-length filters | none | D25 | Identity ≥ 97% removes 93% of R9.4.1 reads and every read of the two R9.4 published mocks; an aligned-length minimum of 1,000-1,400 bp changed nothing useful (§9) |
| Unknown-taxon component | none | D22, D10 | It removed reads of a species missing from the database at 0.9-1.8 times the rate of in-database reads, and reported 14-35% of reads as unknown on mocks whose species are all in the database (§9; `docs/design/unknown_component.md`) |
| Estimator | VB (`em` available) | D2 | VB and EM agree to three decimals on the simulations (§4.3) and within 0.002 in group L1, 0-1 false-positive group apart, on the real samples (§10) |
| Prior, convergence | α₀ = 1, tolerance 1e-7, ≤ 1,000 iterations | F05 | The estimator runs to its fixed point; α₀ = 1 pseudo-read in total matters only for species with a few reads (section 5). Not varied in the benchmark |
| Pruning and floor | Emu's 1/R, 10/R rule; `--min_abundance` 0.001 | D26, D31 | Pruning off reported 94-409 species on real mocks; floors of 1e-4 or 0 at 50k reads raised recall by 0.05-0.07 with 2-8 more false-positive species (§4.3, §5, §6) |
| Intervals | off (`--credible_intervals`) | D7, D27 | 95% coverage 93% / 81% / 38% at 1k / 10k / 50k reads (section 6, §8) |
| `shared_read_fraction` | always | D27 | AUC 0.75-0.84 for flagging misses (§8) |
| Reference database | exact GTDB r226 | D23 | Clustered variants failed through cross-genus identical sequences: species L1 0.52 (99.5% clusters) and 1.24 (ssu_reps) against 0.13 (§4.2) |
| GTDB names | suffixes kept (`--merge-gtdb-suffixes` to merge) | D4 | GTDB treats suffixed species as distinct. Cost on GTDB: siblings split (*L. monocytogenes* / *L. monocytogenes_B* in 70 of 93 samples), +0.4 to +2.3 false-positive species on even simulations scored on merged names; real samples within 0.003 in L1 (§13) |
| Species names | whole (`--binomial_names` to cut) | D4, F17 | Cutting to two words merged distinct species with longer names (`Candidatus` names, NCBI subspecies) |
| Identical-16S column | always | D24 | Species with identical 16S are not identifiable: *S. argenteus* estimated at 0 against a realised 50% of its pair, *E. fergusonii*'s reads given to *E. coli*; species L1 0.131 against group L1 0.010 on the simulations (§4.5, §10) |

**Phase 7b confirmation** (§13). `nextflow run` of the released code with every
default, on the 13 real samples and 80 simulated samples at 1,000 and 10,000 reads,
with the database built with `--merge-gtdb-suffixes`: the same species as the Phase 7a
reference on 93 of 93 samples, species L1 within 0.00045 (real) and 0.00043
(simulations), assigned reads within 3 per 10,000, and every 1k and 10k figure of §4.1
unchanged at three decimals.

## 9. Reference database

`bin/build_gtdb_db.py` builds the database from GTDB's `ssu_all_<release>.fna.gz`
(`--cmd build-database` builds one from any FASTA, sequence-to-taxon table and NCBI or
custom taxonomy):

- **Sequences.** Every 16S sequence of at least 900 bp (`--min-length`) is kept, with
  exact duplicates removed within each species only (MD5 of the uppercase sequence).
  The database keeps every distinct copy of each species, so intragenomic variants and
  strains are all represented.
- **Species.** One tax ID per species name, numbered in alphabetical order, so a
  rebuild gives the same IDs. Lineages come from GTDB's headers.
- **GTDB names are kept with their suffixes by default** (D4, F16): `Listeria
  monocytogenes_B` and `Pseudomonas_E` stay distinct from `Listeria monocytogenes` and
  `Pseudomonas`, as GTDB treats them. `--merge-gtdb-suffixes` strips the suffixes,
  merges those species and genera into one tax ID each, and lists the merges in
  `merged_species.tsv`.
- **GTDB r226.** With `--merge-gtdb-suffixes`: 232,447 sequences over 59,037 species with
  a sequence (81,173 species in `taxonomy.tsv`); 1,529 merged names from 4,787 GTDB
  species, 26 across genera. With the default: 233,611 sequences over 62,295 species with
  a sequence (84,596 in `taxonomy.tsv`).
- **The database is very uneven** (§2, merged build): median 1 sequence per species,
  95th percentile 8, 99th 36, maximum 14,661 (*Escherichia coli*); 39,384 species have
  one sequence and the top 1% of species hold 45.7% of the sequences. This redundancy
  protects the exact database: a true species' many copies fill the -N slots (§4.2).

**Identical-16S groups** (D24; `identical_16s.tsv`). Species of one genus that share at
least one identical sequence are linked, and linked species form a group (genera are
compared without GTDB's suffixes). A sequence that a species of another genus also
carries still links the species of each genus that carry it. On r226: 4,335 species in
1,264 groups in the merged build, 5,215 species in 1,513 groups with suffixes kept
(where suffixed siblings are separate species). Reads from such species score the same
against each other, so the model cannot split them (section 10). The pipeline copies
each species' partners into the `identical_16s_species` column.

**Cross-genus identical sequences** (D23; `cross_genus_identical.tsv`). 1,360 r226
sequences occur in species of more than one genus (1,288 across families), carried by
2,080 species of the merged build. They are almost certainly contaminant contigs in
GTDB assemblies: GTDB's only 16S for *Nonomuraea* sp040654455 is a contig 100% identical
over 1,550 bp to an *S. aureus* 16S (§2). They stay in the database: on the exact
database at -N 3 they carried no false-positive abundance, while in the clustered
variants, which lose the true species' redundancy, they took up to 82% of the estimate
(§4.2). Removing them is future work: removing all of them leaves 355 species with no
sequence (D23).

## 10. Limitations

- **16S copy number is not corrected.** Abundances are shares of 16S reads. Species
  with more 16S copies per genome are over-represented relative to their cell counts.
- **Reads from taxa missing from the database go to the nearest species.** No filter or
  unknown-taxon component ships (D22, D25): in clade-exclusion tests, 58-66% of a
  missing species' reads went to species absent from the sample, or mostly to its
  sibling when that was present (49-57%), and neither identity or aligned-length filters
  nor an unknown component separated them from in-database reads (§9).
- **Species with identical 16S cannot be separated.** Reads that match two species
  equally are split by abundance, so the more abundant species takes them (*S. argenteus*
  at 0 against 50% in its pair; *E. coli* 0.77 of the *E. coli* / *E. fergusonii* pair
  against 0.50 at -N 3, R10.4.1; §4.5). The `identical_16s_species` column names such
  species; read their total as the group's abundance.
- **Near-identical species are only partly separable**, and less so on older chemistry:
  pairs that differ split close to the truth on R10.4.1 reads, while on R9.4.1 11-19% of
  their reads went to the wrong member, and Zymo's *Bacillus* spread over seven species
  of the *B. subtilis* complex on R9.4 mocks (§4.5, §6). These species are not in
  `identical_16s.tsv`; `shared_read_fraction` flags them.
- **Detection is set by the pruning thresholds**: 11-20 reads at 10,000 reads, 51-100 at
  50,000 (section 4). Lowering `--min_abundance` recovers rare species at the cost of
  false positives.
- **The intervals are sampling intervals.** Their coverage falls to 38% at 50,000 reads
  (section 6).
- **GTDB's suffixed siblings split** with the default build (section 8). Use
  `--merge-gtdb-suffixes` for names that follow standard nomenclature.
- **The benchmark's limits** (§12): one simulated community, real data subsampled to
  10,000 reads, and species (not genus) exclusion. The gut mock's reads deviate from its
  theoretical composition in a pattern of amplification bias, which the read-level
  model does not address (§6).

## 11. What the evidence supports

| Claim | Supported? |
|---|---|
| Reproducible Nextflow workflow with pinned containers, tests and CI | Yes |
| GTDB-based reference database builder | Yes, for the exact database; its species definition and cross-genus identical sequences are documented above |
| phyloseq-ready consolidated outputs | Yes |
| A Bayesian (Dirichlet) formulation of Emu's model | Yes, as a description; Emu is the source of the read score and the EM baseline |
| More accurate than Emu | Only as "NanoVI with its exact GTDB database" against Emu with its default database; on the same database they agree (§10) |
| Uncertainty quantification | As sampling intervals plus an ambiguity flag, not as calibrated intervals (§8) |
| Species-level resolution of close relatives | Not for species with identical 16S; partly for near-identical species on R10.4.1 (§4.5) |

Per-claim detail: [benchmarks/m10_claims.md](benchmarks/m10_claims.md).
