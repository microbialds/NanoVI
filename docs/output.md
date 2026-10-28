# NanoVI output

What each command writes, column by column. The examples come from the e2e fixture
(`nextflow run main.nf -profile test,docker`, 630 reads: 600 from six species in the
database, 30 that are not 16S). How the numbers are computed is in
[methods.md](methods.md).

## `--cmd abundance`

```
<output_dir>/
├── sample_outputs/<sample>_rel-abundance.tsv
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
only the two placeholder rows, with all its reads in `Unassigned`; the run log warns, naming it.

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
| `index_cache/md5_cache.tsv` | The MD5 of the database FASTA (and of `taxonomy.tsv` when the identical-16S table is built), with each file's path, size and modification time in milliseconds: a later launch reads the MD5 from here instead of the file while those match. Deleting it only makes the next launch read the files again |
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

## `bin/build_gtdb_db.py`

Builds a database directory from a GTDB SSU FASTA (`ssu_all_<release>.fna.gz`); see the
README for the options. It writes:

| File | Content |
|---|---|
| `species_taxid.fasta` | GTDB 16S sequences of at least `--min-length` bp (900), identical sequences within a species kept once, headers `>taxid:db_name:n`. Tax IDs number the species names in alphabetical order, so a rebuild from the same file gives the same IDs |
| `taxonomy.tsv` | `tax_id`, `species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom` (GTDB's domain), one row per species in the input, including species with no sequence that passes the length filter. Names keep GTDB's clade suffixes unless `--merge-gtdb-suffixes` is given |
| `identical_16s.tsv` | One row per species that shares at least one identical sequence with another species of its genus (genera compared with GTDB's suffixes stripped, so `Pseudomonas_E` counts as `Pseudomonas`). Columns: `tax_id`, `species`, `genus`; `group` (`G1`, `G2`, …, largest first: the species of the genus linked through shared sequences, directly or through a chain); `group_size`; `identical_16s_species` (the species it shares a sequence with directly, `; `-separated, the value the pipeline copies into its column). A sequence that a species of another genus also carries still links the species of each genus that carry it; only the link across genera is left out |
| `cross_genus_identical.tsv` | One row per sequence found in species of more than one genus. Columns: `md5` (of the uppercase sequence), `species_count`, `sequence_ids` and `tax_ids` of the carriers, and the `species`, `genera`, `families` and `phyla` they span (`;`-separated). In GTDB these are almost always contaminant contigs in an assembly. They stay in the database: on the exact database at `-N 3` they carried no false-positive abundance in the benchmark |
| `merged_species.tsv` | Only with `--merge-gtdb-suffixes`: one row per merged name that joins two or more GTDB species with a sequence in the database. Columns: `tax_id`, `species` (the merged name), `gtdb_species_merged` (how many), `gtdb_species` and `gtdb_genera` (suffixes kept, `;`-separated), `ranks_in_conflict` (ranks at which their lineages still differ after stripping, `,`-separated), `lineage_from` (the GTDB species, first in the input, whose lineage `taxonomy.tsv` keeps) |

GTDB r226 (`ssu_all_r226`): with `--merge-gtdb-suffixes`, 232,447 sequences over 59,037
species with a sequence (81,173 species in `taxonomy.tsv`), and `merged_species.tsv`
lists 1,529 names built from 4,787 GTDB species, 26 of them across genera;
`identical_16s.tsv` lists 4,335 species in 1,264 groups. With the default (suffixes
kept), 233,611 sequences over 62,295 species with a sequence (84,596 species in
`taxonomy.tsv`), and `identical_16s.tsv` lists 5,215 species in 1,513 groups: GTDB
siblings such as *Listeria monocytogenes* and *L. monocytogenes_B* are separate species
of one genus there. `cross_genus_identical.tsv` lists the same 1,360 sequences in both
builds. Example rows (default build):

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
