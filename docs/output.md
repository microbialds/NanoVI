# NanoVI output

What each command writes, column by column. The examples come from the e2e fixture
(`nextflow run main.nf -profile test,docker`, 630 reads: 600 from six species in the
database, 30 that are not 16S). How the numbers are computed is in
[methods.md](methods.md).

## `--cmd abundance`

```
<output_dir>/
├── sample_outputs/<sample>_rel-abundance.tsv
├── sample_outputs/<sample>_read_assignments.tsv.gz    (--read_assignments)
├── sample_outputs/<sample>_groups.tsv                 (--ambiguity_groups)
├── sample_outputs/<sample>_cell-abundance.tsv         (--cell_abundance)
├── consolidated_output/
│   ├── otu_table_abundance.tsv
│   ├── otu_table_counts.tsv        (not written when no sample has counts)
│   ├── tax_table.tsv
│   └── sample_metadata.tsv
├── qc/fastplong/<sample>_report.html.gz, <sample>_report.json.gz
├── filtered/<sample>_filtered.fastq.gz        (--keep_files)
├── bam/<sample>_vi_alignments.bam             (--keep_files)
├── index_cache/                               (unless --index_cache_dir points elsewhere)
└── pipeline_info/timeline.html, report.html, trace.txt, software_versions.yml
```

### Per-sample table: `sample_outputs/<sample>_rel-abundance.tsv`

Tab-separated, one row per species the sample is estimated to contain, then two fixed
rows, `Below detection` and `Unassigned`, always last and always present so that tables
from different samples have the same placeholder rows. Species rows are sorted by name. A sample
with no read aligned to the database (a blank, or every read shorter than `--min_length`) has
only the two placeholder rows: `Unassigned` holds the reads that passed the length filter (0 when
none did); the run log warns, naming it.

| Column | When | Meaning |
|---|---|---|
| `species` | always | Species name as the database's `taxonomy.tsv` gives it: whole, GTDB clade suffixes included (`Listeria monocytogenes_B`). With `--binomial_names`, cut to genus and epithet (`Candidatus` names keep three words) and rows that then share a name summed. A tax ID that `taxonomy.tsv` lacks, or that has no species name there, is reported as `Unknown taxid <id>` with a warning |
| `abundance` | always | Posterior mean relative abundance of the species among the reads with at least one scored alignment, renormalised over the reported species, so the species rows sum to 1. It is a share of 16S reads, not of cells (no copy-number correction). `0` on the two placeholder rows |
| `genus` … `superkingdom` | always | Lineage from `taxonomy.tsv` (`genus`, `family`, `order`, `class`, `phylum`, `superkingdom`). On the placeholder rows and on `Unknown taxid` rows with blank ranks, each rank repeats the row's label, so collapsing to a higher rank keeps them apart |
| `estimated_counts` | `--keep_counts` (default on) | Reads attributed to the species: `abundance` × (reads with a scored alignment − reads below detection), rounded by largest remainder so the species rows add up exactly. `Below detection` holds the reads whose candidate species were all pruned; `Unassigned` the reads with no scored alignment (no alignment, or only alignments using an error operation absent from the sample's error profile). All rows together add up to the reads that passed the length filter |
| `identical_16s_species` | always | The species of the same genus that share at least one identical 16S sequence with this species in the database, `; `-separated; empty for most species and on the placeholder rows. 16S cannot separate these species, so the split between them is not evidence for either (see below) |
| `shared_read_fraction` | always | The share of this species' expected reads whose responsibility it splits with another reported species (another species takes at least 1% of the read). 0 means every read is the species' own; a high value means the abundance rests on reads another species explains as well. Empty on the placeholder rows |
| `ci_lower`, `ci_upper` | `--credible_intervals` (VB only) | 2.5% and 97.5% quantiles of the species' mean-field posterior marginal: 95% **sampling** intervals, not calibrated uncertainty (they covered the true abundance in 93%, 81% and 38% of cases at 1,000, 10,000 and 50,000 reads). Empty on the placeholder rows |

On the e2e fixture (with `--credible_intervals`; lineage columns left out):

| species | abundance | estimated_counts | identical_16s_species | shared_read_fraction | ci_lower | ci_upper |
|---|---|---|---|---|---|---|
| Escherichia coli | 0.263519 | 158 | | 0.351 | 0.229 | 0.299 |
| Pseudomonas aeruginosa | 0.100111 | 60 | | 0.0 | 0.077 | 0.125 |
| Shigella flexneri | 0.036536 | 22 | | 0.771 | 0.023 | 0.053 |
| Staphylococcus aureus | 0.399612 | 240 | | 0.0 | 0.361 | 0.439 |
| Streptococcus mitis | 0.050141 | 30 | | 0.067 | 0.034 | 0.069 |
| Streptococcus pneumoniae | 0.150081 | 90 | | 0.101 | 0.123 | 0.180 |
| Below detection | 0.0 | 0 | | | | |
| Unassigned | 0.0 | 30 | | | | |

The truth is 0.25 *E. coli* and 0.05 *S. flexneri*. 51 of their 180 reads score the
same against both species, so the model gives most of them to the more abundant
*E. coli*, and *S. flexneri*'s `shared_read_fraction` of 0.77 says its estimate rests on
such reads. No species in this small database shares an identical 16S with another, so
`identical_16s_species` is empty; on GTDB r226, for example, *Staphylococcus argenteus*
lists *Staphylococcus aureus; Staphylococcus schweitzeri; Staphylococcus singaporensis*.

**Where `identical_16s_species` comes from.** From the database's `identical_16s.tsv`
(below). For a database directory without one, the pipeline's IDENTICAL_16S step builds
it once from `species_taxid.fasta` and `taxonomy.tsv` and keeps it in the index cache as
`species_taxid.<fasta md5>.<taxonomy md5>.identical_16s-<script md5>.tsv` (first 12 hex
digits of each MD5). An `identical_16s.tsv` in `--db` is used as it is when the run's taxonomy is
`<db>/taxonomy.tsv`; with another `--taxonomy_tsv` it is ignored and the table is built
from the FASTA and that taxonomy. With `--db_index`
and no FASTA in `--db`, the column is left empty. Species whose 16S genes are nearly,
but not exactly, identical (for example the *Bacillus subtilis* complex on R9.4.1
reads) are not listed.

**Summed rows.** With `--binomial_names` (or a name repeated in a custom taxonomy), a
row sums several tax IDs: `abundance` and `estimated_counts` add up,
`identical_16s_species` is the union (without the row's own name), `shared_read_fraction`
counts only the reads split with species of other rows (a read split between the summed
tax IDs is the row's own), and the interval uses their summed Dirichlet parameter.

### Per-read assignments: `sample_outputs/<sample>_read_assignments.tsv.gz`

Only with `--read_assignments`. A gzipped, tab-separated long table with one line per
read and row of the per-sample table. It lists every read that passed the length filter:
the scored reads in the order of the alignment file, then the reads with no scored
alignment.

| Column | Meaning |
|---|---|
| `read_id` | The read's name in the FASTQ |
| `species` | A row of the per-sample table, named as there (`--binomial_names` merges rows here too), or `Below detection` or `Unassigned` |
| `probability` | The share of the read that goes to the row: its responsibility in the final fit (φ after the last pruning round for `vb`, the EM posterior for `em`), summed over the row's tax IDs, at full precision |

A line is written when its probability is at least 0.01, the same share that
`shared_read_fraction` uses (there a read is shared when the other rows together take
0.01 of it, so a read split 0.985 / 0.005 / 0.005 / 0.005 is shared but has one line
here; and there the split is judged among the reported rows only, so a read's
`Below detection` share does not make it shared). The lines of one read therefore sum to 1 less what they leave out, which is under
0.01 for each other row. Within a read, rows are listed by decreasing probability and
`Below detection` comes last. Two rows of the table with the same species name and
different lineages (possible in a custom taxonomy) appear under that one name.

- `Below detection` is the read's share on species that are not reported, 1 − Σ of its
  responsibility on reported species. It is exactly 1 for the table's `Below detection`
  reads, all of whose species were pruned. With `vb` it is otherwise 0, unless pruning
  stopped at its iteration limit (the log warns) and `--min_abundance` then dropped
  species. With `em` above 10,000 reads the share can be a fraction: EM refits on the
  species at or above `--min_abundance`, while its responsibilities, as for
  `shared_read_fraction`, come from the fit before.
- `Unassigned` (probability 1) marks the reads with no scored alignment, as counted in the
  table's `Unassigned` row.
- Summed over the file, a species' probabilities give its expected reads, less what falls
  below 0.01. The table's `estimated_counts` are the abundance × assigned reads, rounded.

e2e fixture (`--read_assignments true`): 630 reads in 713 lines. 547 reads have one
line, and 83 split between two rows, 72 of them between *E. coli* and *S. flexneri*.
The first read below leaves out its 0.009 on *S. flexneri*:

```
read_id                    species                   probability
3_Escherichia_coli_0       Escherichia coli          0.9910712754011362
8_Shigella_flexneri_29     Shigella flexneri         0.7659940793128952
8_Shigella_flexneri_29     Escherichia coli          0.2340059206871047
non16S_15                  Unassigned                1.0
```

### Ambiguity groups: `sample_outputs/<sample>_groups.tsv`

Only with `--ambiguity_groups`. The per-sample table's species rows grouped where the
sample cannot separate them, so that a pair of species whose reads are split by
abundance is also reported by its sum, the quantity the reads support. The species table
does not change. Two rows a and b are linked when they are in the same family of the
database's `taxonomy.tsv` and the reads they share carry at least 0.25 of either one's
expected reads, max(f_a(b), f_b(a)) ≥ 0.25, where f_a(b) is the share of a's expected
reads on reads that b takes at least 0.01 of (responsibilities from the final fit, `vb`
or `em`, summed over each row's tax IDs, so `--binomial_names` merges rows first). A row
whose family is blank, or whose tax ID is missing from the taxonomy, links to nothing; a
taxonomy without a family column links no rows (RUN_VI warns). Families are compared with
the ranks above them, so a name used in two orders, such as `Family_XI`, is two families.
The family rule keeps out pairs such as *Akkermansia* and *Faecalibacterium* that share
reads only through a database sequence found in more than one genus. A group is a set of
rows connected by links, so every species row is in exactly one group and every group
lies inside one family; a row that links to nothing is a group of one. f is counted per
partner, so it is not `shared_read_fraction`, which counts a read as shared when all the
other rows together take 0.01 of it.

| Column | When | Meaning |
|---|---|---|
| `group` | always | `G1`, `G2`, … by decreasing abundance; then `Below detection` and `Unassigned` as in the species table |
| `members` | always | The group's species rows, named as in the species table, `; `-separated, by decreasing abundance |
| `n_members` | always | Number of species rows (empty on the placeholder rows) |
| `abundance` | always | Sum of the members' abundances |
| `lowest_common_rank`, `lowest_common_taxon` | always | The lowest rank at which the members' lineages agree and the taxon there (`species` and the name for a group of one); never above the family, since a group lies inside one family |
| `estimated_counts` | `--keep_counts` | Sum of the members' `estimated_counts`; the placeholder rows as in the species table, so the column adds up to the reads that passed the length filter |
| `identical_16s_relatives` | always | Species of the members' genera that share an identical 16S with a member in the database (`identical_16s_species`) and are not reported in the sample, `; `-separated: species the group may also contain, whose reads would then have gone to a member. The fit cannot form groups with them, since they are not reported: *Staphylococcus argenteus*, for example, was pruned beside *S. aureus* in every simulated community that had both at 5% each. Empty when no member has one, or without an identical-16S table |
| `max_shared` | always | The largest max(f_a(b), f_b(a)) between two members: how strongly the group is held together. Empty for a group of one |
| `shared_read_fraction` | always | As in the species table, judged between groups: the share of the group's expected reads whose responsibility is split with rows outside the group. For a group of one it equals the species row's value |
| `ci_lower`, `ci_upper` | `--credible_intervals` (VB only) | From the group's summed Dirichlet parameter, Beta(Σ α, A − Σ α), as for a summed species row. For a group of one they equal the species row's |
| `cell_abundance`, `copy_number`, `copy_number_level` | `--cell_abundance` | The sum of the members' cell-level shares (below), the cell-weighted mean of their copy numbers and the coarsest level any member's came from, over the members that have one. Empty for a group none of whose members has a copy number, and on the placeholder rows |

The 0.25 and the family rule come from the 2026 benchmark ([NanoVI-paper `ambiguity_2026.md`](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/ambiguity_2026.md)): on simulated communities the group table's L1 against the truth summed over the same groups was 0.075, against the species table's 0.132 (0.010 when the listed relatives count towards their group), and 0.326 against 0.352 on 13 real mock communities. The cost: about 1.5 of the 12 species that the species table had within 25% of the truth went into a group per simulated sample. Without the family rule, 7 of the 18 groups on the real mocks joined a false-positive species to a mock species of another phylum or class; the rule removes them, with L1 within 0.001, and leaves those false positives as rows of their own (3.5 per real mock against 2.8 without the rule and 5.2 in the species table).

The group table is per sample: groups are formed from each sample's own reads, so two
samples can group the same species differently, and `combine-outputs`, `collapse-taxonomy`
and the consolidated tables work on the species tables only. The per-read file
(`--read_assignments`) keeps its species rows; a group's reads are its members' lines.

e2e fixture (`--ambiguity_groups true`; lineage and `identical_16s_relatives` columns
left out, the latter is empty since the fixture's 11 sequences are all distinct):

| group | members | n_members | abundance | estimated_counts | max_shared | shared_read_fraction |
|---|---|---|---|---|---|---|
| G1 | Staphylococcus aureus | 1 | 0.399612 | 240 | | 0.0 |
| G2 | Escherichia coli; Shigella flexneri | 2 | 0.300055 | 180 | 0.771 | 0.0 |
| G3 | Streptococcus pneumoniae | 1 | 0.150081 | 90 | | 0.101 |
| G4 | Pseudomonas aeruginosa | 1 | 0.100111 | 60 | | 0.0 |
| G5 | Streptococcus mitis | 1 | 0.050141 | 30 | | 0.067 |
| Below detection | Below detection | | 0.0 | 0 | | |
| Unassigned | Unassigned | | 0.0 | 30 | | |

*S. flexneri*'s estimate rests on reads it shares with *E. coli* (f = 0.771: 16.8 of its
21.8 expected reads), so the two are one group, at 0.300 against a truth of 0.30 (0.25 and
0.05); the group's lowest common rank is the family Enterobacteriaceae. *S. pneumoniae* and
*S. mitis*, which the fixture splits correctly, share 0.101 of *S. pneumoniae*'s reads and
stay apart.

### Cell-level abundance: `sample_outputs/<sample>_cell-abundance.tsv`

Only with `--cell_abundance`, which needs `copy_numbers.tsv` in `--db` and the database's
own `taxonomy.tsv` (a database built without copy numbers, or another `--taxonomy_tsv`,
stops the run before any task starts; `--copy_numbers` and `--rrndb_pantaxa` are
build-database options and are ignored here, with a warning). The species table gives each species'
share of 16S reads, which follows its share of 16S genes; this table divides each species
row's abundance by its 16S copy number and renormalises, so that a species with ten copies
per genome no longer weighs ten times one with a single copy. The species table does not
change.

The file starts with two comment lines (`#`; read it with `comment="#"` in pandas, or
`comment.char = "#"` in R): what the table is, and what it leaves out. The `Below
detection` and `Unassigned` reads have no species and so no copy number, and a species row
none of whose tax IDs has a copy number in `copy_numbers.tsv` has none either; they are
left out, and the second line gives their reads and the read-level abundance of those
species. The shares are therefore of the cells of the species listed.

| Column | Meaning |
|---|---|
| `species` | As in the species table, in the same order |
| `abundance` | Cell-level share: the row's abundance / its copy number, renormalised over the rows listed (sums to 1) |
| `read_abundance` | The row's abundance in the species table (all its tax IDs; a tax ID without a copy number in a row that has one is named in the second comment line) |
| `copy_number` | 16S copies per genome. A row that sums several tax IDs (`--binomial_names`, or a name repeated in a custom taxonomy) takes the cell-weighted mean of their copy numbers, Σ abundance / Σ (abundance / copy number), so its cell share is the sum of theirs; its tax IDs without a copy number are left out of it (and counted in the second comment line) |
| `copy_number_level` | Where the copy number comes from: `species` (the species' own value in `copy_numbers.tsv`), or the `genus` or `family` median it was filled from; for a summed row, the coarsest level among its tax IDs |
| `genus` … `superkingdom` | Lineage, as in the species table |

Copy numbers from a genus or family median carry an error the species value does not
(in GTDB r226, filling a species' known value from the others' median is off by 0.45
copies on average at genus level, 0.81 at family level, and within one copy for 91% and
82% of species). Strains of a species also differ: *E. coli*'s complete genomes mostly
have 7 copies, but some 1 or 17.

e2e fixture (`--cell_abundance true`, copy numbers from
`tests/test_data/e2e/copy_number_sources.tsv`; lineage left out):

```
# Cell-level abundance: read-level abundance / 16S copy number (the database's copy_numbers.tsv), renormalised over the rows listed; ...
# Not included: Below detection (0 reads) and Unassigned (30 reads), and 0 species without a copy number (read-level abundance 0).
```

| species | abundance | read_abundance | copy_number | copy_number_level |
|---|---|---|---|---|
| Escherichia coli | 0.203986 | 0.263519 | 7.0 | species |
| Pseudomonas aeruginosa | 0.135615 | 0.100111 | 4.0 | species |
| Shigella flexneri | 0.028282 | 0.036536 | 7.0 | species |
| Staphylococcus aureus | 0.360888 | 0.399612 | 6.0 | species |
| Streptococcus mitis | 0.067923 | 0.050141 | 4.0 | species |
| Streptococcus pneumoniae | 0.203306 | 0.150081 | 4.0 | species |

The fixture's reads are drawn in fixed read-level proportions, so the cell shares they
imply are those proportions over the copy numbers: *S. aureus* 0.361, *S. pneumoniae*
0.203, *E. coli* 0.194, *P. aeruginosa* 0.136, *S. mitis* 0.068, *S. flexneri* 0.039.
With `--ambiguity_groups` the *E. coli* + *S. flexneri* group gets 0.232 (the sum of its
members) and copy number 7.

The table is per sample, like the group table: the consolidated tables,
`collapse-taxonomy` and `combine-outputs` use the species tables.

In the 2026 benchmark ([NanoVI-paper `copy_number_2026.md`](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/copy_number_2026.md)) the cell table's species L1 against the cell-level
composition was 0.130 on 40 simulated communities whose reads were drawn in proportion
to cells × copies (the read-level table: 0.343; 0.116 with each simulated genome's own
copy number), and 0.367 on 12 Zymo mocks against the cell numbers of their manuals'
Table 1 (read-level table 0.438). On the mocks this only confirms the copy numbers: Zymo
derives its cell shares from the same 16S shares and copy numbers. Replacing the
community species' own values with their genus medians raised L1 by 0.106 and 0.033.

### Consolidated tables: `consolidated_output/`

Built from every sample's table, in the layout phyloseq imports (`otu_table`,
`tax_table`, `sample_data`). Rows are species, sorted by name, then `Below detection`
(when any sample has it) and `Unassigned`; columns are samples. Rows that share a species
name are summed; the interval, `identical_16s_species` and `shared_read_fraction`
columns are left out.

| File | Content |
|---|---|
| `otu_table_abundance.tsv` | `species` and one column per sample: `abundance`, 0 where the sample lacks the species |
| `otu_table_counts.tsv` | The same with `estimated_counts` (the placeholder rows hold their reads). Written only when at least one sample has counts; samples without them get 0, with a warning |
| `tax_table.tsv` | `species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom`. A species whose lineage differs between samples gets its most common lineage. The placeholder rows have blank ranks |
| `sample_metadata.tsv` | `sample_id`; `total_abundance` (sum of the species rows, 1.0); `n_taxa_detected` (species with abundance above 0); `assigned_counts` (reads with a scored alignment, below-detection reads included); `unassigned_counts`. The two counts are 0 when the sample has no counts |

e2e fixture: `sample_metadata.tsv` reads `e2e_mock  1.0  6  600  30`.

### Other abundance outputs

| Path | Content |
|---|---|
| `qc/fastplong/<sample>_report.html.gz`, `.json.gz` | fastplong's reports (gzipped): reads, bases, length and quality distributions before and after filtering. `summary.after_filtering.total_reads` in the JSON is the number of reads the pipeline goes on with |
| `filtered/<sample>_filtered.fastq.gz` | With `--keep_files`: the reads that passed the length filter (and `--min_quality`, when set) |
| `bam/<sample>_vi_alignments.bam` | With `--keep_files`: minimap2's alignments, unsorted: per read the primary, up to `--N` secondary and any supplementary alignments, unmapped reads included |
| `index_cache/<fasta stem>.<md5>.<type>.k<k>.minimap2-<version>.mmi` | The minimap2 index built with `-x <type> -k <kmer_size>`, and its `.versions.yml`. A later run with the same database, settings and minimap2 version reuses it. Not built with `--db_index`. A `-stub-run` adds `.stub` to the name |
| `index_cache/...identical_16s-<md5>.tsv` | The identical-16S table built for a database without one, and its `.versions.yml` (above) |
| `index_cache/md5_cache.tsv` | The MD5 of the database FASTA (and of `taxonomy.tsv` when the identical-16S table is built), with each file's path, size, modification time in milliseconds and, for a local file, inode and status-change time: a later launch reads the MD5 from here instead of the file while those match (a file replaced in place, even with its size and date kept, gets a new inode or status-change time). Deleting it only makes the next launch read the files again |
| `pipeline_info/timeline.html`, `report.html`, `trace.txt` | Nextflow's timeline, report and trace; the trace also records each task's requested `cpus`, `memory` and `time` |
| `pipeline_info/software_versions.yml` | Per process, the tool versions it ran (fastplong; minimap2 and samtools; Python and numpy, pandas, scipy, pysam, biopython), plus NanoVI and Nextflow |

## `--cmd build-database`

Writes `<output_dir>/<db_name>/`, the directory to pass as `--db`:

| File | Content |
|---|---|
| `species_taxid.fasta` | One record per distinct sequence of each species (identical sequences and reverse complements within a species kept once), headers `>taxid:db_name:n` |
| `taxonomy.tsv` | With `--ncbi_taxonomy`: `tax_id`, `species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom`, one row per species-level tax ID in `--seq2tax`, sorted by tax ID. With `--taxonomy_list`: that file's header (its first column renamed `tax_id`) and its rows for the tax IDs in `--seq2tax`, in its order, so give it these rank columns and species-level tax IDs |
| `identical_16s.tsv` | As written by `bin/build_gtdb_db.py` (below) |
| `cross_genus_identical.tsv` | As written by `bin/build_gtdb_db.py` (below) |
| `cross_genus_removed.tsv`, `cross_genus_emptied.tsv` | With a `--cross_genus_rule` other than `keep`: as written by `bin/build_gtdb_db.py` (below). For `most-copies`, a sequence's copies are the input sequences its entry stands for (identical sequences and reverse complements of a species in `--sequences`) |
| `copy_numbers.tsv` | With `--copy_numbers` (a TSV with `species` and `copy_number` columns, matched on the species names of `taxonomy.tsv`; other columns are ignored) or `--rrndb_pantaxa` (rrnDB's pan-taxa statistics for NCBI taxonomy, the `median` of the species row whose `taxid` is the tax ID; needs `--ncbi_taxonomy`, whose tax IDs are NCBI's): the columns below (`bin/build_gtdb_db.py`). `n_genomes` is empty for a table and rrnDB's `childcount` for rrnDB. Without either, no file, and `--cell_abundance` cannot use the database |

## `bin/build_gtdb_db.py`

Builds a database directory from a GTDB SSU FASTA (`ssu_all_<release>.fna.gz`); see the
README for the options. It writes:

| File | Content |
|---|---|
| `species_taxid.fasta` | GTDB 16S sequences of at least `--min-length` bp (900), identical sequences within a species kept once, headers `>taxid:db_name:n`. Tax IDs number the species names in alphabetical order, so a rebuild from the same file gives the same IDs |
| `taxonomy.tsv` | `tax_id`, `species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom` (GTDB's domain), one row per species in the input, including species with no sequence that passes the length filter. Names keep GTDB's clade suffixes unless `--merge-gtdb-suffixes` is given |
| `identical_16s.tsv` | One row per species that shares at least one identical sequence with another species of its genus (genera compared with GTDB's suffixes stripped, so `Pseudomonas_E` counts as `Pseudomonas`). Columns: `tax_id`, `species`, `genus`; `group` (`G1`, `G2`, …, largest first: the species of the genus linked through shared sequences, directly or through a chain); `group_size`; `identical_16s_species` (the species it shares a sequence with directly, `; `-separated, the value the pipeline copies into its column). A sequence that a species of another genus also carries still links the species of each genus that carry it; only the link across genera is left out |
| `cross_genus_identical.tsv` | One row per sequence found in species of more than one genus. Columns: `md5` (of the uppercase sequence), `species_count`, `sequence_ids` and `tax_ids` of the carriers, and the `species`, `genera`, `families` and `phyla` they span (`;`-separated). In GTDB most are contaminant contigs in an assembly, but real strains carry some of them (below). With the default `--cross-genus-rule keep` they stay in the database: on the exact database at `-N 3` they carried no false-positive abundance in the benchmark. With another rule the file still lists the sequences as built, since they are what the rule acted on |
| `cross_genus_removed.tsv` | Only with a `--cross-genus-rule` other than `keep`: one row per entry the rule removed from `species_taxid.fasta`. Columns: `sequence_id`, `tax_id`, `species`, `genus`, `md5` (as in `cross_genus_identical.tsv`), `rule`, `kept_in` (the genus, as in `taxonomy.tsv` with GTDB's suffixes stripped, in which the sequence was kept; empty when the rule removed every copy). `remove` removes every copy; `best-genus` keeps the sequence only in the genus with the most entries carrying it; `most-copies` keeps it only in the genus whose entries stand for the most input sequences of at least `--min-length` (about one per genome), the genus with the most entries in the database breaking a tie. A tie the rule cannot break removes every copy. Genera are compared with GTDB's suffixes stripped |
| `cross_genus_emptied.tsv` | Only with a `--cross-genus-rule` other than `keep`: one row per species left with no sequence. Columns: `tax_id`, `species`, `genus`, `sequences_removed`. The species stay in `taxonomy.tsv`, as species with no sequence of `--min-length` do, but no read can be assigned to them |
| `copy_numbers.tsv` | Only with `--gtdb-metadata` (or `--copy-numbers`, a table as for `--cmd build-database`): one row per species of `taxonomy.tsv`. Columns: `tax_id`, `species`; `copy_number` (16S copies per genome, empty when none); `level` (`species`: the median `ssu_count` over the species' `Complete Genome` assemblies with at least one 16S gene; `genus` or `family`: the median of the species values of its genus, or else of its family, ranks compared with the ranks above them; `none`); `n_genomes` (the complete genomes behind a species value); `source` (`gtdb_complete_genomes:` and the metadata files, or `table:` and the table). `ssu_count` is the number of 16S genes GTDB finds in the genome (HMMER with Barrnap's models), the genome's records in `ssu_all`; draft assemblies are not used because their rRNA operons collapse (in r226, 72% of drafts count 0 or 1; complete genomes 5.3 on average). A blank genus or family fills nothing. `--cell_abundance` reads it. A rebuild without a source removes an earlier `copy_numbers.tsv` |
| `merged_species.tsv` | Only with `--merge-gtdb-suffixes`: one row per merged name that joins two or more GTDB species with a sequence in the database. Columns: `tax_id`, `species` (the merged name), `gtdb_species_merged` (how many), `gtdb_species` and `gtdb_genera` (suffixes kept, `;`-separated), `ranks_in_conflict` (ranks at which their lineages still differ after stripping, `,`-separated), `lineage_from` (the GTDB species, first in the input, whose lineage `taxonomy.tsv` keeps) |

GTDB r226 (`ssu_all_r226`): with `--merge-gtdb-suffixes`, 232,447 sequences over 59,037
species with a sequence (81,173 species in `taxonomy.tsv`), and `merged_species.tsv`
lists 1,529 names built from 4,787 GTDB species, 26 of them across genera;
`identical_16s.tsv` lists 4,335 species in 1,264 groups. With the default (suffixes
kept), 233,611 sequences over 62,295 species with a sequence (84,596 species in
`taxonomy.tsv`), and `identical_16s.tsv` lists 5,215 species in 1,513 groups: GTDB
siblings such as *Listeria monocytogenes* and *L. monocytogenes_B* are separate species
of one genus there. `cross_genus_identical.tsv` lists the same 1,360 sequences in both
builds.

With `--gtdb-metadata` (`bac120_metadata_r226.tsv.gz`, `ar53_metadata_r226.tsv.gz`),
`copy_numbers.tsv` gives the default build's 84,596 species: 11,971 (14.2%) their own
value, 31,049 a genus median, 21,629 a family median, 19,947 none (merged build: 11,124,
29,664, 20,469, 19,916 of 81,173). The species values agree with rrnDB 5.10 (99.1% of the
43,167 genomes in both are equal) and give Zymo's copy numbers for 19 of the 20 species of
its D6300 and D6331 mock communities (*Bifidobacterium adolescentis*: 4, where Zymo's
Table 1 implies 5). Abundance-weighted, the species with their own value held 98.6% of the
assigned abundance of 13 real mocks and 59% of 117 stool samples (25% genus, 15% family,
0.8% none).

With a `--cross-genus-rule`, `identical_16s.tsv` is built from the cleaned sequences, so
it loses the links that a removed sequence made; `taxonomy.tsv` and the tax IDs do not
change. GTDB r226, default build (with `--merge-gtdb-suffixes` in brackets): `remove`
takes out 4,281 (4,147) entries of 2,256 (2,080) species and leaves 413 (355) species
with no sequence; `best-genus` 3,266 (3,286) entries, 224 (207) species left empty;
`most-copies` 2,210 (2,203) entries, 138 (133) species left empty. Most of these
sequences are carried by one species in each genus, so `best-genus` ties on 1,095
(1,117) of the 1,360 sequences and removes every copy of them; `most-copies` removes
every copy of 4.

The rule stays `keep` by default (D38). In the 2026 benchmark
([NanoVI-paper `cross_genus_2026.md`](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/cross_genus_2026.md)) each
rule reports 0.1-0.4 fewer false-positive species per sample, but makes species L1 on
the real mock communities worse on 12-13 of 13 samples (+0.004 to +0.006): the mock
strains carry some of these sequences (*E. coli* loses 73-133 entries), so they are not
all contaminants. `most-copies` costs least. Example rows (default build):

```
identical_16s.tsv
70420  Staphylococcus argenteus  Staphylococcus  G79  8  Staphylococcus aureus; Staphylococcus schweitzeri; Staphylococcus singaporensis

merged_species.tsv (--merge-gtdb-suffixes)
42456  Listeria monocytogenes  6  Listeria monocytogenes;Listeria monocytogenes_B;…;Listeria_A monocytogenes_A  Listeria;Listeria_A    Listeria monocytogenes
```

`bin/identical_16s.py --fasta species_taxid.fasta --taxonomy taxonomy.tsv --output
identical_16s.tsv [--cross-genus cross_genus_identical.tsv]` writes the two reports for an
existing database.

## `--cmd collapse-taxonomy`

`<output_dir>/<input name>-<rank>.tsv` (for `S1_rel-abundance.tsv` at genus,
`S1_rel-abundance-genus.tsv`): the ranks from `--rank` up, with `abundance` (and
`estimated_counts`, when the input has it) summed over the rows that share them. The
placeholder rows stay separate, since their ranks carry their label. The other columns
are left out.

## `--cmd combine-outputs`

`<output_dir>/vi-combined-<rank>.tsv`: one row per taxon at `--rank` (with the ranks
above it), one `abundance` column per `<sample>_rel-abundance.tsv` in `--input_dir`
(named after the sample, in file-name order), 0 where a sample lacks the taxon, sorted by
the rank column. A table without the rank column stops the command with an error that
names the file, and a `--input_dir` with no `*_rel-abundance.tsv` stops it before it runs.
