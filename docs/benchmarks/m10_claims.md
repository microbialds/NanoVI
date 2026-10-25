# What the Phase 7a data supports, by M10 claim

One page. Numbers, sources and caveats are in `methods_2026.md` (sections cited as §).
Real = 13 mocks at 10,000 reads (ONT 2025.09 sup, Emu paper, preprint); simulated = the
20-species community at 1k-50k reads.

| # | Claim | Verdict | What the data shows |
|---|---|---|---|
| 1 | Reproducible Nextflow workflow with pinned containers, tests and CI | Supported (Phases 1-6) | Not a benchmark question. Every benchmark arm runs the pipeline's own steps and images; versions and parameters are recorded (§1). |
| 2 | GTDB-based reference database builder | Supported for the exact database only | The exact GTDB r226 database is what makes NanoVI more accurate than Emu (claim 5). Its clustered variants fail: c995 and reps add 6-25 false-positive species per real sample, through GTDB 16S sequences that are identical across genera (likely contamination; on the simulations reps puts 37-82% of the estimate on them) (§2, §4.2). The exact database is protected only by its redundancy, and less so at -N 50. Document the species definition and the cross-genus identical sequences (M6). |
| 3 | phyloseq-ready consolidated outputs | Supported (Phase 2) | Not a benchmark question. |
| 4 | A Bayesian (Dirichlet) formulation of Emu's model | Supported as a description | VB and Emu's EM give the same abundances everywhere: group L1 within 0.002 on every real and simulated sample, 0-1 false-positive group apart (§4.3, §10). The formulation's value is the posterior (claim 6), not better point estimates. Credit Emu. |
| 5 | More accurate than Emu | Supported only as "NanoVI with its GTDB database" | As shipped, NanoVI has lower error than Emu on all 13 real mocks (identical-16S group L1 0.333 vs 0.415; species L1 0.352 vs 0.568) and on the simulations (species L1 lower by 0.09; group L1 equal). On the same database, both ways, the two give the same L1 (within 0.013 on Emu's database) and Emu reports 1-2 fewer false-positive groups per sample (§10). The gain comes from the database, not the model. Kraken2 + Bracken is far worse (group L1 +0.48 real; 21 more false-positive groups). <!-- pending: wf-16s, NanoCLUST; matched-SAM check of the false-positive gap --> |
| 6 | Uncertainty quantification | Not supported as calibrated intervals; supportable as sampling intervals plus a warning | 95% intervals cover the truth in 93% / 81% / 38% of cases at 1k / 10k / 50k reads (species outside the pairs) and 42% for pair members at 10k. Mean-field, read bootstrap and an exact collapsed Gibbs sampler agree, so the intervals are not too narrow for the model: they miss because of systematic read-assignment error that does not shrink with depth (§8). `shared_read_fraction` flags the misses (AUC 0.75-0.84). |
| 7 | Species-level resolution of close relatives | Not supported for species with identical 16S | Pairs whose members share identical 16S copies are not identifiable: *S. argenteus* is estimated at 0 and *E. fergusonii*'s reads go to *E. coli* (§4.5). Pairs that differ split close to truth on R10.4.1, with 11-19% of reads moved on R9.4.1. Zymo's *Bacillus* spreads over 7 species of the *B. subtilis* complex on R9.4 (§6). A species missing from the database is absorbed by its relatives; neither the unknown component nor identity or length filters detect it (§9). Detection needs 11-20 reads at 10k and 51-100 at 50k, set by the thresholds, not the model (§5). |

## Suggested framing

NanoVI extends Emu's read likelihood with a Bayesian abundance model, an exact GTDB
r226 reference and a reproducible Nextflow workflow. With that reference it is more
accurate than Emu with Emu's default database. With the same reference the two agree,
and the remaining difference is in false positives. It reports sampling intervals and
a per-species ambiguity signal (`shared_read_fraction`), not calibrated uncertainty.
Species whose 16S genes are identical are reported, or should be reported, as groups.

## Also needed for the manuscript (W3-W6, W9)

- W3, runtime: <!-- pending: timing.tsv --> §11.
- W4, same-database comparisons both ways: §10 (Emu on NanoVI's database; NanoVI on
  Emu's database).
- W5, false positives per tool: §10, `comparators_per_sample.tsv`.
- W6, detection limit: §5.
- W9, source data: every figure in `figures/` has its data in `data/` under the same name.
