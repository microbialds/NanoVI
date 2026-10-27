# Manuscript figures

`make_figures.py` draws the manuscript figures from the benchmark tables in
`docs/benchmarks/data/` and nothing else: every plotted value is read from those files
(columns described in `docs/benchmarks/methods_2026.md`). The only constants drawn are
reference marks (the nominal 0.95 coverage, zero error, the y = x diagonal). `style.py`
holds the shared style (decision D12).

```bash
# all figures
python docs/paper/make_figures.py --out <output folder>
# one figure (full name or prefix)
python docs/paper/make_figures.py --out <output folder> --only fig3
```

Needs matplotlib (>= 3.10, for `hatchcolor`), pandas and numpy; no seaborn. The manuscript figures were drawn with the benchmark environment's matplotlib 3.11.2, pandas 3.0.6 and numpy 2.4.6 (Python 3.11).

For each figure the script writes `<name>.pdf` and `<name>.svg` (vector, text kept as
text: TrueType/Type 42 fonts in the PDF, `svg.fonttype = none`), `<name>.png` (600 dpi)
and `<name>_data.tsv`: the plotted rows and columns, with a `source_file` column (and a
`panel` column for multi-panel figures).

## Style (D12)

- Width 85 mm (single column) or 175 mm (double column); height as needed, at most 200 mm.
- Arial if installed, otherwise DejaVu Sans; 7 pt text, 8 pt bold lower-case panel labels
  at the top-left of each panel; 0.6 pt axes; no top or right spines; no titles.
- Okabe-Ito palette. Tools and databases keep one colour in every figure: NanoVI on GTDB
  r226 #0072B2, NanoVI on Emu's database #56B4E9, Emu on its database #D55E00, Emu on GTDB
  #E69F00, wf-16s minimap2 #009E73, wf-16s kraken2 #CC79A7, Kraken2 + Bracken #000000,
  NanoCLUST #F0E442 (dark edge).
- Figures with no tool comparison (1, 3-6) show NanoVI on GTDB r226 only, so their other
  encodings reuse colours: in figures 1, 4 and 5 the NanoVI blue for the points, and read
  depth (1k, 10k, 50k) as a light-to-dark ramp of that blue; figure 3 gives its eight
  expected species the eight Okabe-Ito colours; figure 6 colours the filter families.
- Simulated samples (figures 1b and 4a): circle for even, triangle for log-staggered
  communities; filled for R10.4.1, open for R9.4.1.
- Species names are in italics in tick labels and legends.

## Figures

| Figure | Width | Source files (docs/benchmarks/data/) |
|---|---|---|
| `fig1_database` | double | `db_audit_ssu_all_seqs_per_species_hist.tsv`, `db_audit_ssu_all_summary.tsv`, `fig_db_N.tsv` |
| `fig2_comparators` | double | `fig_comparators.tsv` |
| `fig3_zymo_composition` | double | `real_species.tsv` |
| `fig4_pairs_lod` | double | `pairs_ties.tsv`, `fig_lod.tsv` |
| `fig5_intervals` | double | `fig_intervals.tsv`, `intervals_coverage_by_shared_read_fraction.tsv`, `intervals_warning_auc.tsv` |
| `fig6_clade` | single | `fig_clade.tsv` |
| `fig7_runtime` | double | `fig_timing.tsv` |

**fig1_database.** (a) GTDB r226 species by number of database sequences, in log2 bins
(1, 2-3, 4-7, ...; log axes), with the median, 95th percentile and maximum read from
`db_audit_ssu_all_summary.tsv`. Only species with at least one sequence are in the
histogram. (b) Species L1 (left) and genus L1 (right) for the arms exact -N 3/10/50,
c995 -N 3 and reps -N 3: one point per cell sample (four: even and log-staggered, R10.4.1
and R9.4.1, 10k reads, seed 1) and the mean as a bar. `fig_db_N.tsv` is the reference
setting of the cell grid (primary profile, Emu score, VB, floor 1e-3, pruning on); its
means are not the table in methods §4.2, which averages `ablation_cell_db_N.tsv` over the
other cell layers.

**fig2_comparators.** Real mocks (`kind == real`): (a) identical-16S group L1 (`l1_group`)
and (b) false-positive groups per sample (`group_fp_count`, symlog axis) per tool and
database, one point per sample (jittered vertically; circles for Zymo community
standards, squares for gut mocks) and the mean as a bar. Tools are ordered by mean group
L1, best at the top. All tools ran on the 13 real mocks except Emu on GTDB r226 (5) and
NanoCLUST (the 12 Zymo mocks); those labels say how many. NanoVI's EM (`nanovi_em`) is
left out: it matches VB on every sample (methods §10).

**fig3_zymo_composition.** Stacked bars of the expected Zymo 16S composition (one bar: the
`expected` values of the `in_truth` rows, identical for every Zymo sample; the script
stops if they differ) and of the NanoVI estimate (`arm == reference`) for SRR23636353,
SRR10391201 and each of the six MSPlus sup replicates (one bar each, not their mean).
Segments: *B. spizizenii* (the Zymo strain in GTDB) first, then every *Bacillus* species
not in the truth in greys (alternate ones hatched; *B. subtilis* is one of them), the other
seven expected species, the four MSPlus ATCC spike-ins as one hatched white segment
(`in_truth` but no expected share: presence only), and every other non-truth species as
"other false-positive species". The estimate bars are the raw estimate (they sum to 1 and
include the spike-ins); nothing is renormalised.

**fig4_pairs_lod.** (a) For each near-neighbour pair at -N 3 and 10k reads, the VB share of
the first-named member within the pair minus its realised share, reads_a / (reads_a +
reads_b), one point per seed for each chemistry and shape, the mean over seeds as a bar.
(b) Detection rate (`detected`) against the species' realised reads (`reads_bin`, in file
order) for arm `exact_N3_primary_emu_vb_f0.001`, species outside the pairs; one line per
depth (1k, 10k, 50k), solid for R10.4.1 and dashed (open markers) for R9.4.1, the two
chemistries offset slightly so both stay visible where they coincide.

**fig5_intervals.** (a) Coverage of the 95% intervals (`coverage_design`) against depth for
the mean-field, collapsed Gibbs and read-bootstrap intervals; solid for species outside
the pairs, dashed (open markers) for pair members; the methods are offset slightly on x
because they nearly coincide. (b) Mean-field coverage by `shared_read_fraction` bin at
each depth, with the AUC of `shared_read_fraction` as a predictor of misses
(`intervals_warning_auc.tsv`, mean-field rows) in the legend. Dotted line: nominal 0.95.

**fig6_clade.** Share of the excluded species' reads removed against the share of
in-database species' reads removed, for every clade-exclusion arm: colour by family (no
filter, identity, aligned length, unknown component), circle for R10.4.1 and open square
for R9.4.1, arms of one family joined in threshold order; dashed y = x. Both axes are
symlog (linear below 0.001) because several arms remove nothing. The strictest arm of
each filter is labelled.

**fig7_runtime.** CPU seconds (`cpu_s`, medians of 3 repeats) for the two timing samples.
For each -N and k, NanoVI's bar is align + sam_to_bam + nanovi_likelihood + nanovi_vb, and
Emu's is the same alignment + `emu` on the same SAM; Emu as shipped (`mode emu_shipped`) is
one undivided bar. (a) Totals, alignment in grey and the post-alignment part in the tool's
colour (totals printed at the bar ends). (b) Post-alignment steps only, on a linear axis.
NanoVI's EM (`nanovi_em`) and the index build (`mode index`, reported apart in methods
§11) are not drawn.
