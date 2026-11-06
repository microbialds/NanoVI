# NanoVI workflow

How the processes of each `--cmd` connect and which parameter switches each step or output
on. The one-figure overview is [`images/nanovi_overview.svg`](images/nanovi_overview.svg);
the parameters are in the [README](../README.md#parameters), the files in
[output.md](output.md) and the model in [methods.md](methods.md).

**How to read the diagrams.** Every box with an upper-case name is a process or subworkflow
named as in `main.nf` and `subworkflows/local/`; rounded boxes are files. An edge label is
either the condition under which the edge exists, written as `main.nf`, the subworkflow or
the module tests it (`params.ambiguity_groups`), or, for an edge that always exists, what it
carries. Boxes with a solid green border run with the default parameters, two of them
only when the run needs them: INDEX_DB is skipped with `--db_index` or when the index
cache already holds the index, and IDENTICAL_16S runs only when the run has no
identical-16S table to use (so not for a database from `bin/build_gtdb_db.py` or
`--cmd build-database`, which ships one, unless `--taxonomy_tsv` is another file). Boxes
with a dashed border need a parameter.

`tests/check_workflow_diagram.sh` (run in CI) fails when a process or subworkflow that
`main.nf` or a subworkflow includes is not named in a diagram or table on this page, or when
the page shows one that the code no longer includes.

---

## `--cmd abundance` (the default)

`main.nf` checks the inputs, reads the samples (a `.csv` samplesheet checked against
`assets/schema_input.json`, or a folder of FASTQ files) and calls `ABUNDANCE`; `REPORT`
runs after `pipeline_info/software_versions.yml` is written, because it reads it.

```mermaid
flowchart TD
    Reads(["FASTQ per sample<br/>samplesheet .csv or folder"]) -->|"params.input"| ABUNDANCE
    ABUNDANCE["ABUNDANCE<br/>subworkflows/local/abundance.nf"] -->|"each sample"| FILTER_READS

    FILTER_READS["FILTER_READS<br/>fastplong: length 500-2,000 bp<br/>mean quality with --min_quality"] -->|"filtered reads"| ALIGN_READS
    FILTER_READS -->|"params.min_quality != null"| Quality(["--mean_qual filter"])

    Db[("--db directory<br/>species_taxid.fasta, taxonomy.tsv")] -->|"!params.db_index"| INDEX_DB["INDEX_DB<br/>minimap2 -x map-ont -k 21<br/>storeDir: index cache"]
    DbIndex(["--db_index .mmi"]) -->|"params.db_index"| ALIGN_READS
    INDEX_DB -->|"index"| ALIGN_READS
    Db -->|"no identical_16s.tsv in --db,<br/>or --taxonomy_tsv is not the database's"| IDENTICAL_16S["IDENTICAL_16S<br/>bin/identical_16s.py<br/>storeDir: index cache"]
    Db -->|"identical_16s.tsv, used as is"| WRITE_OUTPUT
    IDENTICAL_16S -->|"identical-16S table"| WRITE_OUTPUT

    ALIGN_READS["ALIGN_READS<br/>minimap2 -ax map-ont -N 3 -p 0.9<br/>into samtools view -b"] -->|"BAM"| LIKELIHOOD
    ALIGN_READS -->|"params.keep_files"| Bam(["bam/"])
    FILTER_READS -->|"params.keep_files"| Filtered(["filtered/"])
    LIKELIHOOD["LIKELIHOOD<br/>bin/read_likelihood.py<br/>error profile, Emu's read score"] -->|"read_scores.npz"| RUN_VI
    LIKELIHOOD -->|"params.read_assignments:<br/>read names kept in the npz"| RUN_VI
    LIKELIHOOD -->|"read_counts.json"| WRITE_OUTPUT

    Tax(["taxonomy.tsv<br/>--taxonomy_tsv, else in --db"]) -->|"rows of the table"| RUN_VI
    Tax -->|"lineages"| WRITE_OUTPUT
    Tax -->|"when IDENTICAL_16S runs"| IDENTICAL_16S
    RUN_VI["RUN_VI<br/>bin/estimate_abundance.py<br/>params.estimator: vb or em"] -->|"abundance JSON"| WRITE_OUTPUT
    RUN_VI -->|"params.read_assignments"| ReadAssign(["sample_outputs/#lt;sample#gt;_read_assignments.tsv.gz"])
    RUN_VI -->|"params.ambiguity_groups:<br/>groups in the JSON"| WRITE_OUTPUT

    CopyNumbers(["copy_numbers.tsv in --db"]) -->|"params.cell_abundance"| WRITE_OUTPUT
    WRITE_OUTPUT["WRITE_OUTPUT<br/>bin/write_abundance_table.py"] -->|"always; ci_lower, ci_upper with<br/>params.credible_intervals"| SpeciesTable(["sample_outputs/#lt;sample#gt;_rel-abundance.tsv"])
    WRITE_OUTPUT -->|"params.ambiguity_groups"| Groups(["sample_outputs/#lt;sample#gt;_groups.tsv"])
    WRITE_OUTPUT -->|"params.cell_abundance"| Cell(["sample_outputs/#lt;sample#gt;_cell-abundance.tsv"])
    SpeciesTable -->|"every sample, collected"| GENERATE_PHYLOSEQ_TABLES
    GENERATE_PHYLOSEQ_TABLES["GENERATE_PHYLOSEQ_TABLES<br/>bin/generate_phyloseq_tables.py"] -->|"otu_table_counts.tsv with<br/>params.keep_counts"| Consolidated(["consolidated_output/"])

    Versions(["versions.yml of every task"]) -->|"collectFile"| SoftwareVersions(["pipeline_info/software_versions.yml"])
    SoftwareVersions -->|"params.cmd == 'abundance' && params.report"| REPORT["REPORT<br/>bin/nanovi_report.py"]
    FILTER_READS -->|"params.report: fastplong JSON"| REPORT
    SpeciesTable & Groups & Cell & ReadAssign & Consolidated -.->|"params.report"| REPORT
    Metadata(["--report_metadata"]) -->|"params.report_metadata"| REPORT
    REPORT -->|"params.report"| ReportDir(["report/: index.html, nanovi_report.pdf, figures/"])

    classDef on fill:#E8F4EC,stroke:#009E73,stroke-width:2px,color:#111
    classDef opt fill:#F7F7F7,stroke:#777,stroke-width:1.5px,stroke-dasharray:5 4,color:#111
    classDef data fill:#FFF8E6,stroke:#E69F00,stroke-width:1.5px,color:#111
    class ABUNDANCE,FILTER_READS,INDEX_DB,IDENTICAL_16S,ALIGN_READS,LIKELIHOOD,RUN_VI,WRITE_OUTPUT,GENERATE_PHYLOSEQ_TABLES on
    class REPORT,Quality,Bam,Filtered,ReadAssign,Groups,Cell,ReportDir,Metadata,DbIndex opt
    class Reads,Db,Tax,CopyNumbers,SpeciesTable,Consolidated,Versions,SoftwareVersions data
```

`INDEX_DB` and `IDENTICAL_16S` keep their outputs with `storeDir` in `--index_cache_dir`
(default `<output_dir>/index_cache`) under names that record the database, so a later run
that finds the file skips the task. `IDENTICAL_16S` runs only when the run has no table to
use: `--db` has no `identical_16s.tsv`, or `--taxonomy_tsv` is not `<db>/taxonomy.tsv`;
without a FASTA (`--db_index` alone) the `identical_16s_species` column is left empty.
`--cell_abundance` stops the run at launch when `--db` has no `copy_numbers.tsv` or when
`--taxonomy_tsv` is not the database's own.

---

## `--cmd build-database`

```mermaid
flowchart LR
    Seqs(["--sequences FASTA"]) -->|"params.sequences"| BUILD_DATABASE
    Seq2tax(["--seq2tax table"]) -->|"params.seq2tax"| BUILD_DATABASE
    TaxList(["--taxonomy_list TSV"]) -->|"params.taxonomy_list"| BUILD_DATABASE
    Ncbi(["--ncbi_taxonomy dump"]) -->|"params.ncbi_taxonomy"| BUILD_DATABASE
    CopyTable(["--copy_numbers table"]) -->|"params.copy_numbers"| BUILD_DATABASE
    Rrndb(["--rrndb_pantaxa file"]) -->|"params.rrndb_pantaxa<br/>(needs params.ncbi_taxonomy)"| BUILD_DATABASE
    BUILD_DATABASE["BUILD_DATABASE<br/>subworkflows/local/build_database.nf"] -->|"always"| BUILD_DB["BUILD_DB<br/>bin/main.py build-database"]
    BUILD_DB -->|"always"| DbDir(["#lt;db_name#gt;/: species_taxid.fasta, taxonomy.tsv,<br/>identical_16s.tsv, cross_genus_identical.tsv"])
    BUILD_DB -->|"params.cross_genus_rule != 'keep'"| Removed(["cross_genus_removed.tsv,<br/>cross_genus_emptied.tsv"])
    BUILD_DB -->|"params.copy_numbers or params.rrndb_pantaxa"| Cn(["copy_numbers.tsv"])

    classDef cmd fill:#FFFFFF,stroke:#6E7781,stroke-width:1.5px,color:#111
    classDef opt fill:#F7F7F7,stroke:#777,stroke-width:1.5px,stroke-dasharray:5 4,color:#111
    class BUILD_DATABASE,BUILD_DB cmd
    class CopyTable,Rrndb,Removed,Cn opt
```

`main.nf` needs `--sequences`, `--seq2tax` and exactly one of `--taxonomy_list` or
`--ncbi_taxonomy`, and at most one of `--copy_numbers` or `--rrndb_pantaxa`.
`bin/build_gtdb_db.py` builds a GTDB database outside Nextflow with the same files
([output.md](output.md#binbuild_gtdb_dbpy)).

---

## `--cmd collapse-taxonomy`

```mermaid
flowchart LR
    InTsv(["--input_tsv: a per-sample table"]) -->|"params.input_tsv, params.rank"| COLLAPSE_TAXONOMY
    COLLAPSE_TAXONOMY["COLLAPSE_TAXONOMY<br/>subworkflows/local/postprocessing.nf"] -->|"always"| COLLAPSE["COLLAPSE<br/>bin/main.py collapse-taxonomy"]
    COLLAPSE -->|"always"| Collapsed(["#lt;input name#gt;-#lt;rank#gt;.tsv"])

    classDef cmd fill:#FFFFFF,stroke:#6E7781,stroke-width:1.5px,color:#111
    class COLLAPSE_TAXONOMY,COLLAPSE cmd
```

---

## `--cmd combine-outputs`

```mermaid
flowchart LR
    InDir(["--input_dir with *_rel-abundance.tsv"]) -->|"params.input_dir, params.rank"| COMBINE_OUTPUTS
    COMBINE_OUTPUTS["COMBINE_OUTPUTS<br/>subworkflows/local/postprocessing.nf"] -->|"always"| COMBINE["COMBINE<br/>bin/main.py combine-outputs"]
    COMBINE -->|"always"| Combined(["vi-combined-#lt;rank#gt;.tsv"])

    classDef cmd fill:#FFFFFF,stroke:#6E7781,stroke-width:1.5px,color:#111
    class COMBINE_OUTPUTS,COMBINE cmd
```

`main.nf` stops `combine-outputs` before it runs when `--input_dir` holds no
`*_rel-abundance.tsv`.

---

## Parameter gates

Each step or output and the condition under which it runs, as the code tests it.

| Step | Runs when | Where the condition is |
|---|---|---|
| `ABUNDANCE`, `FILTER_READS`, `ALIGN_READS`, `LIKELIHOOD`, `RUN_VI`, `WRITE_OUTPUT`, `GENERATE_PHYLOSEQ_TABLES` | `params.cmd == 'abundance'` (the default) | `main.nf` |
| `FILTER_READS`: mean-quality filter | `params.min_quality != null` | `modules/local/filter_reads.nf` |
| `INDEX_DB` | `!params.db_index` (no prebuilt index) | `subworkflows/local/abundance.nf` |
| `IDENTICAL_16S` | no `identical_16s.tsv` in `--db`, or `--taxonomy_tsv` other than `<db>/taxonomy.tsv`, and a FASTA to build it from | `main.nf`, `subworkflows/local/abundance.nf` |
| `LIKELIHOOD`: read names in the npz | `params.read_assignments` | `modules/local/likelihood.nf` |
| `RUN_VI`: `<sample>_read_assignments.tsv.gz` | `params.read_assignments` | `modules/local/inference.nf` |
| `RUN_VI`: groups in the JSON; `WRITE_OUTPUT`: `<sample>_groups.tsv` | `params.ambiguity_groups` | `modules/local/inference.nf`, `write_output.nf` |
| `WRITE_OUTPUT`: `<sample>_cell-abundance.tsv` | `params.cell_abundance` (and `copy_numbers.tsv` in `--db`) | `main.nf`, `modules/local/write_output.nf` |
| `WRITE_OUTPUT`: `ci_lower`, `ci_upper` | `params.credible_intervals` (needs `--estimator vb`) | `main.nf`, `modules/local/write_output.nf` |
| `WRITE_OUTPUT`: `estimated_counts`; `GENERATE_PHYLOSEQ_TABLES`: `otu_table_counts.tsv` | `params.keep_counts` (default `true`) | `modules/local/write_output.nf` |
| `ALIGN_READS`, `FILTER_READS`: publish `bam/`, `filtered/` | `params.keep_files` | `conf/modules.config` |
| `REPORT` | `params.cmd == 'abundance' && params.report` | `main.nf` |
| `BUILD_DATABASE`, `BUILD_DB` | `params.cmd == 'build-database'` | `main.nf` |
| `BUILD_DB`: `copy_numbers.tsv` | `params.copy_numbers` or `params.rrndb_pantaxa` | `main.nf` |
| `BUILD_DB`: `cross_genus_removed.tsv`, `cross_genus_emptied.tsv` | `params.cross_genus_rule != 'keep'` | `modules/local/build_database.nf` |
| `COLLAPSE_TAXONOMY`, `COLLAPSE` | `params.cmd == 'collapse-taxonomy'` | `main.nf` |
| `COMBINE_OUTPUTS`, `COMBINE` | `params.cmd == 'combine-outputs'` | `main.nf` |

---

## Processes and subworkflows

### Subworkflows

| Subworkflow | Command | Calls |
|---|---|---|
| `ABUNDANCE` | `abundance` | `FILTER_READS`, `INDEX_DB`, `IDENTICAL_16S`, `ALIGN_READS`, `LIKELIHOOD`, `RUN_VI`, `WRITE_OUTPUT`, `GENERATE_PHYLOSEQ_TABLES` |
| `BUILD_DATABASE` | `build-database` | `BUILD_DB` |
| `COLLAPSE_TAXONOMY` | `collapse-taxonomy` | `COLLAPSE` |
| `COMBINE_OUTPUTS` | `combine-outputs` | `COMBINE` |

### Processes

Containers and resources are set per label in `conf/base.config`; with `-profile conda` the
same tools come from bioconda, the Python steps from `containers/environment.yml` and the
report from `containers/report/environment.yml`. `<version>` is NanoVI's version. Every
process also writes a `versions.yml`, collected into `pipeline_info/software_versions.yml`
and not published beside the results.

| Process | Label | Image | Publishes |
|---|---|---|---|
| `FILTER_READS` | `filter_reads` | `quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0` | `qc/fastplong/<sample>_report.{html,json}.gz`; `filtered/<sample>_filtered.fastq.gz` with `--keep_files` |
| `INDEX_DB` | `alignment` | BioContainers mulled minimap2 2.31 + samtools 1.24 (`mulled-v2-66534bcb…`) | the index into `--index_cache_dir` (storeDir) |
| `IDENTICAL_16S` | `vi_database` | `ghcr.io/microbialds/nanovi-python:<version>` | the identical-16S table into `--index_cache_dir` (storeDir) |
| `ALIGN_READS` | `alignment` | as `INDEX_DB` | `bam/<sample>_vi_alignments.bam` with `--keep_files` |
| `LIKELIHOOD` | `vi_stats` | `ghcr.io/microbialds/nanovi-python:<version>` | nothing (the npz and counts JSON stay in the work directory) |
| `RUN_VI` | `vi_inference` | `ghcr.io/microbialds/nanovi-python:<version>` | `sample_outputs/<sample>_read_assignments.tsv.gz` with `--read_assignments` |
| `WRITE_OUTPUT` | `vi_output` | `ghcr.io/microbialds/nanovi-python:<version>` | `sample_outputs/<sample>_rel-abundance.tsv`; `<sample>_groups.tsv` with `--ambiguity_groups`; `<sample>_cell-abundance.tsv` with `--cell_abundance` |
| `GENERATE_PHYLOSEQ_TABLES` | `vi_output` | `ghcr.io/microbialds/nanovi-python:<version>` | `consolidated_output/`: `otu_table_abundance.tsv`, `otu_table_counts.tsv` (with counts), `tax_table.tsv`, `sample_metadata.tsv` |
| `REPORT` | `report` | `ghcr.io/microbialds/nanovi-report:<version>` | `report/` with `--report` |
| `BUILD_DB` | `vi_database` | `ghcr.io/microbialds/nanovi-python:<version>` | `<output_dir>/<db_name>/` |
| `COLLAPSE` | `vi_collapse` | `ghcr.io/microbialds/nanovi-python:<version>` | `<output_dir>/<input name>-<rank>.tsv` |
| `COMBINE` | `vi_combine` | `ghcr.io/microbialds/nanovi-python:<version>` | `<output_dir>/vi-combined-<rank>.tsv` |

Nextflow's own `pipeline_info/timeline.html`, `report.html` and `trace.txt` are written for
every command.
