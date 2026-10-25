# Unknown-taxon component (design note, M5)

Status: prototype behind experimental options (Phase 7a, decision D10). Not used by the pipeline; Phase 7b decides whether to implement it, from the clade-exclusion benchmark in `docs/benchmarks/methods_2026.md`.

## Problem

A read counts as unassigned only when it has no alignment at all. A read from a species missing from the database still aligns to its closest relative, at lower identity, and both estimators (VB and EM) must give it to one of the species it aligns to. It inflates that relative, and when the relative is absent from the sample it creates a false positive. Environmental and clinical samples often contain species with no 16S sequence in any database. Emu has the same limitation; its 3.6 release added hard filters on alignment identity and aligned length (`--min-pid`, `--min-align-len`).

## Model

The mixture gets one more component, *u*, for "a species not in the database":

- z_r ∈ {1, …, S} ∪ {u}, π ~ Dirichlet(α₀/(S+1), …, α₀/(S+1)) over the S + 1 components.
- p(x_r | z_r = s) ∝ exp(ℓ(r, s)) for a database species s, as now (ℓ is Emu's score, or the per-base score of M4).
- p(x_r | z_r = u) ∝ exp(c_r), where c_r is a fixed per-read *floor score*: what the read would score against a reference that differs from it by the sample's typical error rate plus a divergence δ.

So a read whose best alignment scores well above its floor is explained by a database species; a read whose every alignment scores below the floor is better explained by "something else". How much better is not decided read by read: the floor competes through π_u like any other species, so a few poor reads in a sample of good ones still go mostly to their best species, while a sample with many poor reads builds up π_u.

## The floor, per sample

For read r with longest alignment length L_max(r):

- **Emu score**: c_r = L_max(r) × (m + δ) × Σ_op p_op ln p_op, where p_op is the sample's error profile (the shares of I, D, S and X among errors in the primary alignments; the same values the score uses) and m is the sample's median per-read error rate over primary alignments (errors over M + I + D + S columns). This is the score of a virtual alignment of the read's longest alignment length with error rate m + δ and the sample's error mix.
- **Per-base score**: c_r = L_max(r) × [(1 − (m + δ)) ln p_match + (m + δ) Σ_op w_op ln p_op], with the per-base log-probabilities of the score and w_op the error mix.

m adapts the floor to the chemistry (about 6% for R9.4.1, 1% for R10.4.1 SUP) without a parameter; δ is the one parameter. A read from an absent species aligns to its nearest relative with roughly its own errors plus the divergence between the two genes; with δ near the 16S species boundary (about 1.3%, 98.7% identity) the floor catches reads whose nearest database species is a different species, while reads from a database species with ordinary errors stay above it. The benchmark tries δ ∈ {0.01, 0.02, 0.03, 0.05}.

Alternatives not prototyped: a per-read baseline from the read's own quality scores (mean 10^(−Q/10)), which follows read-to-read variation but needs calibrated Q scores and does not see soft clips; and a floor learned per sample (an extra parameter fitted by EM), which the identifiability of π_u and c_r makes fragile.

## Updates

VB (mean field), with E[ln π_k] = ψ(α_k) − ψ(Σ_j α_j):

- φ_rs ∝ exp(ℓ(r, s) + E[ln π_s]) for each candidate s of read r, and φ_ru ∝ exp(c_r + E[ln π_u]); the responsibilities of a read sum to 1 over its candidates and u.
- α_s = α₀/(S+1) + Σ_r φ_rs and α_u = α₀/(S+1) + Σ_r φ_ru.

EM: the same with ln π_k for E[ln π_k], and π_k ∝ Σ_r φ_rk.

In the prototype, *u* is one more column of the sparse likelihood (tax ID 0, `UNKNOWN_TAX_ID`) with one entry per scored read, so both estimators, the pruning loop and the convergence rule apply unchanged. Pruning treats it like any species: when fewer than max(10/R, floor) of the reads go to it, it is pruned and those reads return to their best database species. An implementation could exempt it from pruning.

## Reporting

Prototype: the abundance JSON carries tax ID 0 with its share π_u; the benchmark renormalises the database species over themselves and reports π_u separately.

Proposed for the pipeline (Phase 7b, if adopted): a row `Unknown (not in database)` next to `Below detection` and `Unassigned`, with estimated_counts = π_u × assigned reads (rounded with the other rows) and abundance 0, the database species renormalised over themselves, so that the abundance column keeps its meaning (shares of reads explained by database species). The phyloseq tables would carry it like the other placeholder rows. A sample whose unknown share is large is then visible at a glance.

## Relation to identity and length filters

Emu's filters drop alignments below a fixed identity or aligned length before scoring: a read whose alignments are all dropped becomes unassigned. They are hard cuts at the same threshold in every sample and ignore abundance. The floor is soft (a read below it is shared between u and its species according to the abundances), adapts to each sample's error rate, and keeps the read in the denominator as `Unknown`, not `Unassigned`. The benchmark compares the two on the clade-exclusion sets: where the excluded species' reads go, false positives created, and the effect on species that are in the database.

## Prototype options

- `read_likelihood.py --unknown-delta δ`: writes `unknown_floor` (one value per scored read, in read order) into the scores npz (`likelihood.unknown_floor_scores`).
- `estimate_abundance.py --unknown-component`: adds the candidate for each read and reports its abundance under tax ID 0 (both estimators).
- For comparison: `read_likelihood.py --min-identity` and `--min-aligned-length`, with Emu 3.6's definitions (percent identity 100 × (1 − NM / aligned query bases); aligned query bases exclude clips).

## Benchmark result (Phase 7a)

Source: `docs/benchmarks/methods_2026.md` Section 9, `docs/benchmarks/data/clade_exclusion_summary.tsv`,
`real_arms_paired.tsv`.

- **Novelty detection fails on the clade-exclusion sets.** With one community species
  removed from the database at a time, the component removes the excluded species'
  reads at 1.1-1.8 times the rate at which it removes reads of species in the database
  (δ 0.02: 8% vs 6% on R10.4.1, 31% vs 28% on R9.4.1). Identity and length filters do
  no better. Every excluded species has GTDB relatives within the read error of its
  16S, so per-read evidence cannot separate a missing species from a present relative.
- **On real mocks it acts as a soft filter.** With every species in the database, δ 0.02
  reports 14-35% of the reads as unknown (R10.4.1 sup 14-20%, R9.4 34-35%). The
  database species, renormalised over themselves, keep their species L1 (within 0.01),
  and 1-4 false-positive species per sample disappear.
- **Reporting as proposed above would mislead.** An `Unknown (not in database)` row of
  14-35% on a mock with no novel taxa would be read as novelty. If the component is
  kept for its false-positive trimming, it should be reported as a share of reads
  that fit no database sequence within the sample's error, not as unknown taxa.

Whether to adopt it, and how to name the row, are decisions for Phase 7b.

