# NanoVI methods

How NanoVI turns full-length 16S Nanopore reads into species abundances, what it
takes from Emu and what it changes, the evidence behind each default and option, how the
GTDB reference is built, and what the method cannot do. Numbers come from the 2026
benchmark, [methods_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/methods_2026.md) (cited as §), the three
later reports on the database's cross-genus sequences
([cross_genus_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/cross_genus_2026.md)),
the ambiguity groups ([ambiguity_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/ambiguity_2026.md)) and
the 16S copy-number correction ([copy_number_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/copy_number_2026.md)),
cited by name and section (cross_genus §6), and their source data in
[`docs/benchmarks/data/`](https://github.com/microbialds/NanoVI-paper/tree/main/docs/benchmarks/data), all in the
[NanoVI-paper](https://github.com/microbialds/NanoVI-paper) repository with the
benchmark scripts. Decisions are cited by their roadmap keys
(D2-D42). The output files are described in [output.md](output.md), and every process
with the parameter that switches it on in [workflow.md](workflow.md).

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
output names them, and an optional table reports the species a sample cannot separate
as groups with their summed abundance (section 7). Other options, all off by default,
write each read's species probabilities (section 6), a cell-level table corrected for
16S copy number (section 11) and a run report (section 12), and clean the database of
16S sequences found in more than one genus (section 10); with the defaults none of them
changes the output.

## 1. Pipeline

| Step | Tool | What it does (defaults) |
|---|---|---|
| FILTER_READS | fastplong 0.7.0 | Keeps reads of 500-2,000 bp (`--min_length`, `--max_length`). The quality filter is off (`--disable_quality_filtering`); `--min_quality Q` turns on fastplong's mean-quality filter (`--mean_qual Q`) |
| INDEX_DB | minimap2 2.31 | Indexes the database's `species_taxid.fasta` with `-x map-ont -k 21`, once; the index is cached under a name that records the FASTA's MD5, the preset, k and the minimap2 version |
| IDENTICAL_16S | `bin/identical_16s.py` | Only for a database without `identical_16s.tsv`: the species of each genus that share an identical 16S (section 10), built once and cached beside the index |
| ALIGN_READS | minimap2 2.31, samtools 1.24 | `minimap2 -ax map-ont -N 3 -p 0.9 -K 500M -k 21`: the primary alignment and up to 3 secondary alignments scoring at least 90% of it; BAM output |
| LIKELIHOOD | `bin/read_likelihood.py` | The sample's error profile and each read's score against each candidate species (section 2) |
| RUN_VI | `bin/estimate_abundance.py` | Abundances by variational inference with pruning (sections 3-4); `--estimator em` runs Emu's EM instead (section 5). With `--read_assignments`, each read's probabilities (section 6); with `--ambiguity_groups`, the groups (section 7) |
| WRITE_OUTPUT, GENERATE_PHYLOSEQ_TABLES | `bin/write_abundance_table.py`, `bin/generate_phyloseq_tables.py` | Per-sample tables with lineage, counts, `identical_16s_species` and `shared_read_fraction` (section 6), and the group and cell-level tables when their options are on (sections 7, 11); consolidated phyloseq tables |
| REPORT | `bin/nanovi_report.py` | With `--report`: an HTML site and a PDF from the published tables (section 12) |

[workflow.md](workflow.md) has a diagram of each command, with the parameter that
switches each step or output on, and the container image of each process.

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
matches was tested and not adopted (section 9, D29).

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
as a check on the first iteration, every 100 iterations and on the last (a fall of more
than $10^{-9}|\mathcal{L}|$ between two computed values prints a warning); it does not
decide convergence (F51).

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

## 6. Posterior summaries: `shared_read_fraction`, intervals and per-read probabilities

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

**Per-read probabilities** (`--read_assignments`, off by default; D37). The same final
responsibilities, written per read: for every read that passed the length filter, each
row $g$ of the table with $\sum_{s \in g} \phi_{rs} \ge 0.01$ and that sum, then
`Below detection` $= 1 - \sum_{s \in \mathcal{R}} \phi_{rs}$ (exactly 1 for the reads
whose species were all pruned) and `Unassigned` (probability 1) for the reads with no
scored alignment. The 0.01 is the share `shared_read_fraction` uses, but there a read is
shared when the other rows together take 0.01 of it. Summed over the file, a row's
probabilities give its expected reads, less what falls below 0.01. Nothing else
changes; LIKELIHOOD then also keeps the read names in its scores file. On 200,000
synthetic reads over 2,000 species RUN_VI took 6.8-7.0 s and 199 MB with the file,
against 4.1-4.3 s and 153 MB without (Phase 12, F55).

## 7. Ambiguity groups

Where species have identical or nearly identical 16S genes, their reads score about the
same against each, the fit splits them by abundance, and one member of a pair can take
most of the reads or all of them (section 13). `--ambiguity_groups` (off by default;
D39, F57) reports such rows together, with the sum the reads support, in a separate
per-sample table; the species table does not change.

**Rule.** For two rows $a$ and $b$ of the table, from the final fit (VB's $\phi$ or EM's
responsibilities, summed over each row's tax IDs, so `--binomial_names` merges rows
first),

$$
f_a(b) = \frac{1}{E_a} \sum_r \phi_{ra}\, \mathbf{1}\big[\tilde\phi_{rb} \ge 0.01\big],
\qquad E_a = \sum_r \phi_{ra},
$$

the share of $a$'s expected reads on reads that $b$ takes at least 0.01 of
($\tilde\phi$ as in section 6, among the reported rows). $a$ and $b$ are linked when
they are in the same family of `taxonomy.tsv` (compared with the ranks above it) and
$\max(f_a(b), f_b(a)) \ge \tau = 0.25$; the maximum catches a small row whose estimate
rests on reads it shares with a large one. A group is a connected component of the
links, so every row is in exactly one group (alone when it links to nothing) and every
group lies inside one family. A row with no family links to nothing. A group reports
the summed abundance and counts, its lowest common rank, $\max f$ among its members, its
own `shared_read_fraction` (reads split with rows outside the group) and, with
`--credible_intervals`, an interval from the members' summed $\alpha$. It also lists
the species that share an identical 16S with a member in the database and are not
reported (`identical_16s_relatives`): the fit cannot group them, since pruning removed
them, but the group may contain them.

**Evidence** (ambiguity §3, §4, §7; the same alignments as Phase 13, `-N 3`, merged exact
database). On the 92 simulations the group table's L1, against the truth summed over the
same groups, was 0.075 against the species table's 0.132 (lower on 83, higher on none),
and 0.0096 when the listed relatives' truth counts towards their group; the
*E. coli* / *E. fergusonii*, *Bacillus* and *Streptococcus* groups were within 25% of
their truth in all 152 samples where they formed. The cost: 1.5 of the 12.1 species
per sample that the species table had within 25% of the truth went into a group, and
separable pairs were grouped where one member was much rarer (*Klebsiella* in 22 of 40
log-distributed samples). On the 13 real mocks, L1 fell from 0.352 to 0.326, almost all
on the three R9.4.1 mocks. Without the family rule, 7 of the 18 real-mock groups crossed
a genus, each through a false-positive species that carries one of the cross-genus
sequences of section 10; the rule removes them, changes mean L1 by at most 0.0006 and
leaves those species as rows of their own (3.54 false-positive rows per real mock,
against 2.85 without the rule and 5.15 in the species table). τ = 0.2 scored marginally
lower but hid more right species (1.8 per sample) and formed more cross-genus groups;
τ = 0.3 gave up most of the real-mock gain (0.337). The group intervals covered 0.69 of
the group truths at 10,000 reads, against 0.26 for their members' intervals.

Groups are formed from each sample's reads, so two samples can group the same species
differently; the consolidated tables, `collapse-taxonomy` and `combine-outputs` use the
species tables only.

## 8. What differs from Emu, and why

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
| Extra outputs | | `shared_read_fraction`, `identical_16s_species`, phyloseq tables; optional sampling intervals, per-read probabilities, ambiguity groups, a cell-level table and a run report | D24, D27, D37, D39, D40, D42 |
| Database | rrnDB + NCBI 16S RefSeq, NCBI taxonomy | GTDB r226 `ssu_all`, all distinct sequences per species, GTDB taxonomy (section 10) | The database is what makes NanoVI more accurate than Emu as shipped (§10) |

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

## 9. Defaults and their evidence

Every default below was kept or set after the 2026 benchmark (Phase 7a) and confirmed
on the released code (Phase 7b, §13); the last five rows are the options added in
release 1.2.0, each benchmarked in its own report and off by default.

| Setting | Default | Decision | Evidence |
|---|---|---|---|
| Quality filter | off (length only) | D5, D28 | Q10 and Q15 changed real-mock species L1 by less than 0.008 and the false-positive count by at most one species; on the simulations, thresholds high enough to matter removed reads unevenly across species (R10.4.1 even, genus L1 0.005 → 0.035 at Q33) (§7) |
| Tool versions | fastplong 0.7.0, minimap2 2.31 | D30 | Against fastplong 0.2.2 and minimap2 2.24: the same species on every real sample, species L1 within 0.0004, assigned reads within 3 per 10,000 (§7, §13) |
| Secondary alignments | `-N 3` | D21 | Simulations: species L1 0.130 / 0.143 / 0.182 at -N 3 / 10 / 50 on the exact database (genus L1 0.008 / 0.012 / 0.046), while false-positive species fall (1.6 / 0.5 / 0.25 per sample): at -N 50 GTDB's cross-genus identical copies take alignment slots and move reads between species of the community (§4.2). Real mocks: -N 50 removes 1-3 false-positive species per sample, lowers species L1 on the R9.4 Zymo mocks (−0.05, −0.08), changes it by less than 0.003 on the ONT 2025.09 data, raises it on the Emu gut mock (+0.04) and loses one gut species in five samples (§6). -N 3 is also the cheapest |
| Error profile | primary alignments | D3, D31 | 'All alignments' changed no metric at 1k or 10k reads; species L1 +0.0008 at 50k (§4.3) |
| Read score | Emu's | D29 | A per-base score raised recall on the simulations (10k even: 0.990 vs 0.950) but added false positives on real mocks (L1 −0.004 to +0.012, 0-2 more false-positive species) (§4.4) |
| Identity and aligned-length filters | none | D25 | Identity ≥ 97% removes 93% of R9.4.1 reads and every read of the two R9.4 published mocks; an aligned-length minimum of 1,000-1,400 bp changed nothing useful (§9) |
| Unknown-taxon component | none | D22, D10 | It removed reads of a species missing from the database at 0.9-1.8 times the rate of in-database reads, and reported 14-35% of reads as unknown on mocks whose species are all in the database (§9; [design note](https://github.com/microbialds/NanoVI-paper/blob/main/docs/design/unknown_component.md) in NanoVI-paper) |
| Estimator | VB (`em` available) | D2 | VB and EM agree to three decimals on the simulations (§4.3) and within 0.002 in group L1, 0-1 false-positive group apart, on the real samples (§10) |
| Prior, convergence | α₀ = 1, tolerance 1e-7, ≤ 1,000 iterations | F05 | The estimator runs to its fixed point; α₀ = 1 pseudo-read in total matters only for species with a few reads (section 5). Not varied in the benchmark |
| Pruning and floor | Emu's 1/R, 10/R rule; `--min_abundance` 0.001 | D26, D31 | Pruning off reported 94-409 species on real mocks; floors of 1e-4 or 0 at 50k reads raised recall by 0.05-0.07 with 2-8 more false-positive species (§4.3, §5, §6) |
| Intervals | off (`--credible_intervals`) | D7, D27 | 95% coverage 93% / 81% / 38% at 1k / 10k / 50k reads (section 6, §8) |
| `shared_read_fraction` | always | D27 | AUC 0.75-0.84 for flagging misses (§8) |
| Reference database | exact GTDB r226 | D23 | Clustered variants failed through cross-genus identical sequences: species L1 0.52 (99.5% clusters) and 1.24 (ssu_reps) against 0.13 (§4.2) |
| GTDB names | suffixes kept (`--merge-gtdb-suffixes` to merge) | D4 | GTDB treats suffixed species as distinct. Cost on GTDB: siblings split (*L. monocytogenes* / *L. monocytogenes_B* in 70 of 93 samples), +0.4 to +2.3 false-positive species on even simulations scored on merged names; real samples within 0.003 in L1 (§13) |
| Species names | whole (`--binomial_names` to cut) | D4, F17 | Cutting to two words merged distinct species with longer names (`Candidatus` names, NCBI subspecies) |
| Identical-16S column | always | D24 | Species with identical 16S are not identifiable: *S. argenteus* estimated at 0 against a realised 50% of its pair, *E. fergusonii*'s reads given to *E. coli*; species L1 0.131 against group L1 0.010 on the simulations (§4.5, §10) |
| Cross-genus sequences | kept (`--cross-genus-rule keep`) | D38 | Each of the three rules gave 0.1-0.4 fewer false-positive species per sample, but made species L1 on the real mocks worse on 12-13 of 13 samples (+0.004 to +0.006), because the mock strains carry some of these sequences; the simulations, which exclude every genome with such a sequence, improved by at most 0.002 (cross_genus §5, §6) |
| Per-read file | off (`--read_assignments`) | D37 | An audit file, not an estimate: it changes no other output (section 6) |
| Ambiguity groups | off (`--ambiguity_groups`); τ = 0.25, links within a family | D39 | Group L1 0.075 against species L1 0.132 on the simulations and 0.326 against 0.352 on the real mocks, at the cost of 1.5 right species per simulated sample hidden in a group (section 7; ambiguity §3, §4, §7) |
| Cell-level table | off (`--cell_abundance`); GTDB complete-genome median, then genus, then family median | D40 | Species L1 against the cell composition 0.130 against 0.343 for the read-level table on 40 simulations and 0.367 against 0.438 on 12 Zymo mocks (section 11; copy_number §4, §5) |
| Run report | off (`--report`) | D42 | Reads only what the run publishes; the other outputs do not change (section 12) |

**Phase 7b confirmation** (§13). `nextflow run` of the released code with every
default, on the 13 real samples and 80 simulated samples at 1,000 and 10,000 reads,
with the database built with `--merge-gtdb-suffixes`: the same species as the Phase 7a
reference on 93 of 93 samples, species L1 within 0.00045 (real) and 0.00043
(simulations), assigned reads within 3 per 10,000, and every 1k and 10k figure of §4.1
unchanged at three decimals.

## 10. Reference database

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
against each other, so the model cannot split them (section 13). The pipeline copies
each species' partners into the `identical_16s_species` column.

**Cross-genus identical sequences** (D23, D38; `cross_genus_identical.tsv`). 1,360 r226
sequences occur in species of more than one genus (1,288 across families), carried by
2,080 species of the merged build. Most look like contaminant contigs in GTDB
assemblies: GTDB's only 16S for *Nonomuraea* sp040654455 is a contig 100% identical over
1,550 bp to an *S. aureus* 16S (§2). On the exact database at -N 3 they carried no
false-positive abundance, while in the clustered variants, which lose the true species'
redundancy, they took up to 82% of the estimate (§4.2). `--cross-genus-rule`
(`--cross_genus_rule` for `--cmd build-database`) can clean them: `remove` every copy;
`best-genus`, keep each sequence only in the genus with the most species carrying it; or
`most-copies`, only in the genus whose entries stand for the most input sequences (about
one per GTDB genome); a tie the rule cannot break removes every copy. On the merged
build the three remove 4,147, 3,286 and 2,203 entries and leave 355, 207 and 133
species with no sequence (cross_genus §3). The default stays `keep` (D38): each rule
reported 0.1-0.4 fewer false-positive species per sample, but species L1 on the 13 real
mocks rose from 0.352 to 0.356-0.358, worse on 12-13 of 13 samples, 92-100% of it on mock
species that lost entries (*E. coli* lost 73-133), so real strains carry some of these
sequences; the simulations, which by construction contain none of them, improved by at
most 0.002 (cross_genus §5, §6). With a rule, `identical_16s.tsv` describes the cleaned
database.

**16S copy numbers** (D40; `copy_numbers.tsv`, only with `--gtdb-metadata`). GTDB's
`ssu_count` (16S genes found by HMMER with Barrnap's models, the genome's records in
`ssu_all`) is usable only for complete genomes: draft assemblies collapse their rRNA
operons, and 72% of those in species with complete genomes count 0 or 1. A species'
copy number is the median `ssu_count` over its `Complete Genome` assemblies (14.2% of the
default build's 84,596 species; rrnDB 5.10 matched by accession covers 11.5% and gives
the same value for 99.1% of the 43,167 genomes in both). A species without one takes the
median of the species values of its genus, else of its family (ranks keyed with those
above them): 37% and 26% of species, with a leave-one-out error of 0.45 and 0.81 copies
(91% and 82% within one copy); 24% have none (copy_number §2). The values give Zymo's own
copy numbers for 19 of the 20 species of its D6300 and D6331 mocks. `--cmd build-database`
takes a species table (`--copy_numbers`) or rrnDB's pan-taxa statistics (`--rrndb_pantaxa`)
and fills them the same way.

## 11. Cell-level abundance (16S copy number)

The species table gives each species' share of 16S reads, which follows its share of 16S
genes, not of cells. With `--cell_abundance` (off by default; D40, F58) each sample also
gets a table with, for each row $g$ of the species table with copy number $k_g$,

$$
c_g = rac{\pi_g / k_g}{\sum_{h} \pi_h / k_h} ,
$$

over the rows that have a copy number. A row of several tax IDs (`--binomial_names`), and
with `--ambiguity_groups` a group, takes the cell-weighted mean
$k_g = \sum_{s \in g} \pi_s / \sum_{s \in g} (\pi_s / k_s)$, so its cell share is the
sum of its members'. The `Below detection` and `Unassigned` reads and rows with no copy
number have no cell share; the file's first two lines say how much read-level abundance
that leaves out. The species and consolidated tables do not change.

**Evidence** (copy_number §4-§6). On 40 simulated communities whose reads were drawn in
proportion to cells × copies of one complete genome per species, species L1 against the
cell composition was 0.130 for the cell table against 0.343 for the read-level table
(lower on 40 of 40); with each simulated genome's own copy number it was 0.116, so the
gap between a strain and its species median adds 0.014. Genus L1: 0.027 against 0.265. On
the 12 Zymo mocks with a cell-level composition (their manuals' Table 1) it was 0.367
against 0.438 (lower on 11 of 12); this only confirms the copy numbers, since Zymo derives
its cell shares from the same 16S shares and copy numbers. Replacing the community
species' own values with their genus medians raised L1 by 0.106 on the simulations and
0.033 on the mocks (*P. aeruginosa* 4 filled as 6, *E. faecalis* 4 as 6). In a stool
cohort of 117 samples, 59% of the assigned abundance had a species-level value, 25% a
genus and 15% a family median.

## 12. Run report

`--report` (off by default; D42, F60) adds a last step that reads only what the run
publishes (the species tables, the consolidated tables, fastplong's reports,
`software_versions.yml`, and the group, cell-level and per-read tables when their options
ran) and writes an offline HTML site, a PDF and figures (SVG, PDF, 600 dpi PNG, each with a
CSV of its values): each sample's read fate, the composition at species and genus, the
species flagged by `shared_read_fraction` ≥ 0.5 or a non-empty `identical_16s_species`
(a display rule, not a test), the groups and the cell-level table, a summary of the
per-read file, and the methods with every tool's version and reference. It runs in its
own image, so the Python steps' image and outputs are unchanged, and two renders of the
same results give identical files. The design follows the run report of BugBuster
(gene2dis/BugBuster): one report model in Python, Typst for the PDF, Jinja for the site,
matplotlib for the figures and Vega-Lite for the interactive ones.

## 13. Limitations

- **The species table is read-level.** Abundances are shares of 16S reads, so species
  with more 16S copies per genome are over-represented relative to their cells.
  `--cell_abundance` divides by a copy number (section 11), but most species' value is a
  genus or family median (off by 0.45 and 0.81 copies on average in GTDB r226), strains
  differ from their species' median (*E. coli*'s complete genomes mostly have 7 copies,
  some 1 or up to 19), 24% of GTDB species have no value and are left out, and the
  correction cannot repair species the read-level estimate misses or splits wrongly. The
  mock panel only checks the copy numbers against Zymo's; no mock with measured cell
  counts was tested.
- **Reads from taxa missing from the database go to the nearest species.** No filter or
  unknown-taxon component ships (D22, D25): in clade-exclusion tests, 58-66% of a
  missing species' reads went to species absent from the sample, or mostly to its
  sibling when that was present (49-57%), and neither identity or aligned-length filters
  nor an unknown component separated them from in-database reads (§9). The cross-genus
  rules send more of them to absent species (cross_genus §7).
- **Species with identical 16S cannot be separated.** Reads that match two species
  equally are split by abundance, so the more abundant species takes them (*S. argenteus*
  at 0 against 50% in its pair; *E. coli* 0.77 of the *E. coli* / *E. fergusonii* pair
  against 0.50 at -N 3, R10.4.1; §4.5). The `identical_16s_species` column names such
  species, and `--ambiguity_groups` reports the rows the sample cannot separate with their
  sum (section 7). A group cannot include a species pruned away (*S. argenteus* beside
  *S. aureus*), only list it among `identical_16s_relatives`, and it can hide a species
  that the species table had right (1.5 of 12.1 per simulated sample).
- **Near-identical species are only partly separable**, and less so on older chemistry:
  pairs that differ split close to the truth on R10.4.1 reads, while on R9.4.1 11-19% of
  their reads went to the wrong member, and Zymo's *Bacillus* spread over seven species
  of the *B. subtilis* complex on R9.4 mocks (§4.5, §6). These species are not in
  `identical_16s.tsv`; `shared_read_fraction` flags them, and the group table joins them
  where they share enough reads in the sample.
- **The database keeps 16S sequences found in more than one genus.** Most are contaminant
  contigs, but cleaning them made the real mocks worse (section 10), so they stay by
  default; they can link a false-positive species to a mock species of another family,
  which the family rule keeps out of the groups.
- **Detection is set by the pruning thresholds**: 11-20 reads at 10,000 reads, 51-100 at
  50,000 (section 4). Lowering `--min_abundance` recovers rare species at the cost of
  false positives.
- **The intervals are sampling intervals.** Their coverage falls to 38% at 50,000 reads
  (section 6).
- **GTDB's suffixed siblings split** with the default build (section 9). Use
  `--merge-gtdb-suffixes` for names that follow standard nomenclature.
- **The benchmark's limits** (§12): one simulated community, real data subsampled to
  10,000 reads, and species (not genus) exclusion. The gut mock's reads deviate from its
  theoretical composition in a pattern of amplification bias, which the read-level
  model does not address (§6).

## 14. What the evidence supports

| Claim | Supported? |
|---|---|
| Reproducible Nextflow workflow with pinned containers, tests and CI | Yes |
| GTDB-based reference database builder | Yes, for the exact database; its species definition and cross-genus identical sequences are documented above |
| phyloseq-ready consolidated outputs | Yes |
| A Bayesian (Dirichlet) formulation of Emu's model | Yes, as a description; Emu is the source of the read score and the EM baseline |
| More accurate than Emu | Only as "NanoVI with its exact GTDB database" against Emu with its default database; on the same database they agree (§10) |
| Uncertainty quantification | As sampling intervals plus an ambiguity flag, not as calibrated intervals (§8) |
| Species-level resolution of close relatives | Not for species with identical 16S; partly for near-identical species on R10.4.1 (§4.5) |
| Reporting what 16S cannot separate | Yes, as groups of rows within a family: group L1 0.075 against species L1 0.132 on the simulations, 0.326 against 0.352 on the real mocks, with a cost in species hidden (section 7) |
| Cell-level abundances | As an estimate whose error follows the copy numbers: better than the read-level table against cell compositions (0.130 against 0.343 on simulations), checked on mocks only against Zymo's own copy numbers (section 11) |
| Cleaning the database of cross-genus sequences | Not as a default: fewer false positives, worse real-mock abundances (section 10) |

Per-claim detail: [m10_claims.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/m10_claims.md) (NanoVI-paper).
