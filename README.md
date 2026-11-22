# NanoVI

NanoVI tells you which bacterial and archaeal species are in each of your samples, and in what proportions, from full-length 16S rRNA reads sequenced on Oxford Nanopore. It is a pipeline that you run with one command: [Nextflow](https://www.nextflow.io/) downloads it and runs every step in a container, so there is nothing to compile or install besides Nextflow and Docker (or Singularity, or conda).

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![CI](https://github.com/microbialds/NanoVI/actions/workflows/ci.yml/badge.svg)](https://github.com/microbialds/NanoVI/actions/workflows/ci.yml)

[![What NanoVI does: your FASTQ files of full-length 16S Nanopore reads, one per sample, are filtered to reads of 500 to 2,000 bases (fastplong), matched to the 16S sequences of a GTDB reference database built once (minimap2), and scored to estimate each species' share of the reads (NanoVI's Bayesian model, with Emu's read score); the results are a species table per sample and one for the whole run, ready for R and phyloseq. Optional extras: groups of species that 16S cannot tell apart, abundances corrected for 16S copies per genome, the species each read was assigned to, and an HTML and PDF report.](docs/images/nanovi_overview.svg)](docs/images/nanovi_overview.png)

## At a glance

- **What goes in:** one FASTQ file per sample of full-length 16S rRNA amplicon reads from Oxford Nanopore (for example the 16S Barcoding Kit, basecalled with Dorado), and a 16S reference database that NanoVI builds for you once from [GTDB](https://gtdb.ecogenomic.org/).
- **What comes out:** for each sample, a table of the species found with their relative abundance, estimated read counts and full lineage; for the whole run, the same as a species × samples table with a taxonomy table and sample summary, ready to load into R with phyloseq. On request it also writes groups of species that 16S cannot tell apart, abundances corrected for the number of 16S copies per genome, the species each read was assigned to, and an HTML and PDF report with figures.
- **What you need:** a Linux computer or server (or an HPC cluster) with Java 17 or later, Nextflow 26.04 or later, and Docker, Singularity/Apptainer or conda. You run commands in a terminal; no programming is needed, and you do not need to download this repository.
- **Who it is for:** anyone with 16S Nanopore amplicon data who wants species-level profiles. NanoVI reports abundances as shares of 16S reads (see [what the numbers mean](#what-the-numbers-mean)); species whose 16S genes are identical cannot be told apart by any 16S method, and NanoVI says when that happens.

## Getting started

### 1. Install the requirements

```bash
java -version                     # needs 17 or later; if missing, install a JDK, e.g. from https://adoptium.net
curl -s https://get.nextflow.io | bash
chmod +x nextflow && mkdir -p ~/bin && mv nextflow ~/bin/    # make sure ~/bin is on your PATH
nextflow -version                 # needs 26.04 or later; `nextflow self-update` upgrades it
docker run --rm hello-world       # checks that Docker works for your user
```

Docker: see [docs.docker.com/engine/install](https://docs.docker.com/engine/install/). On Linux, your user must be allowed to run Docker (`sudo usermod -aG docker $USER`, then log out and in again). On a cluster without Docker, use `-profile singularity` instead of `-profile docker` below (see [Profiles](#profiles-docker-singularity-conda-and-clusters)). More on installing Nextflow: [nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html).

### 2. Try NanoVI on the bundled test data

```bash
nextflow run microbialds/NanoVI -r v1.2.0 -profile test,docker --output_dir nanovi_test
```

Nextflow downloads NanoVI (into `~/.nextflow/assets/`) and its container images (about 1.2 GB, once), then runs the whole pipeline on a small test set: 630 simulated reads from six species and an 11-species database. Once the images are downloaded the run takes under a minute. The species table is `nanovi_test/consolidated_output/otu_table_abundance.tsv`; it should report *Staphylococcus aureus* at about 0.40 and *Escherichia coli* at about 0.26 (the expected values are in [tests/test_data/e2e/expected.tsv](tests/test_data/e2e/expected.tsv)).

- `-r v1.2.0` runs that release. Without `-r`, Nextflow runs the newest commit of the `main` branch; give the release so that you can report which version you ran (it is also written to `pipeline_info/software_versions.yml`).
- `microbialds/NanoVI` is the GitHub repository; Nextflow does not mind the capital letters (`microbialds/nanovi` works too).
- `-profile` has a single dash (it is a Nextflow option); NanoVI's own options have two (`--output_dir`).

### 3. Build the reference database (once)

NanoVI's reference database is built from [GTDB](https://gtdb.ecogenomic.org/)'s 16S sequences. GTDB is the database NanoVI is benchmarked with, and the one we recommend; release r226 is the one the benchmarks used. Download its 16S file and, for the optional cell-level abundances, its two metadata files (70 MB, 225 MB and 5 MB), and check them against GTDB's list of checksums:

```bash
mkdir -p gtdb_r226 && cd gtdb_r226
B=https://data.gtdb.ecogenomic.org/releases/release226/226.0
wget $B/genomic_files_all/ssu_all_r226.fna.gz $B/bac120_metadata_r226.tsv.gz $B/ar53_metadata_r226.tsv.gz $B/MD5SUM.txt
grep -E '(ssu_all|bac120_metadata|ar53_metadata)_r226' MD5SUM.txt | sed 's#  .*/#  #' | md5sum -c -   # three "OK" lines
cd ..
```

Then build the database:

```bash
nextflow run microbialds/NanoVI -r v1.2.0 -profile docker \
    --cmd build-database \
    --gtdb_ssu gtdb_r226/ssu_all_r226.fna.gz \
    --gtdb_metadata gtdb_r226/bac120_metadata_r226.tsv.gz,gtdb_r226/ar53_metadata_r226.tsv.gz \
    --db_name db_gtdb_r226 \
    --output_dir databases
```

This takes about 2 minutes and under 1 GB of memory and writes `databases/db_gtdb_r226/` (about 370 MB), the directory you give to every run as `--db`. Leave out `--gtdb_metadata` if you will not use `--cell_abundance`. [Reference databases](#reference-databases) has the options, the files, newer GTDB releases and other databases.

### 4. Describe your samples

Write a samplesheet, a CSV file with a `sample` column (the name you want in the results) and a `fastq` column (the path to that sample's reads):

```csv
sample,fastq
patient01,/data/run1/fastq_pass/barcode01.fastq.gz
patient02,/data/run1/fastq_pass/barcode02.fastq.gz
```

Each sample is one FASTQ file (`.fastq`, `.fq`, optionally gzipped). If your basecaller wrote many small files per barcode, join them first, for example `cat fastq_pass/barcode01/*.fastq.gz > barcode01.fastq.gz`. If all your files are in one folder and their names are good sample names, you can give the folder instead of a samplesheet ([Input](#input-samplesheet-or-folder)).

### 5. Run

```bash
nextflow run microbialds/NanoVI -r v1.2.0 -profile docker \
    --input samplesheet.csv \
    --db databases/db_gtdb_r226 \
    --output_dir results
```

The first run with a database also builds its alignment index (under a minute for GTDB on 8 CPUs) and keeps it in `results/index_cache/`; later runs into the same `--output_dir` reuse it (to share it between projects, see `--index_cache_dir` under [Parameters](#parameters)). If a run stops, add `-resume` to the same command to continue from the finished steps. Running again into the same `--output_dir` replaces the earlier results.

### 6. Read the results

| File | What it is |
|------|------------|
| `results/consolidated_output/otu_table_abundance.tsv` | Species × samples: each species' share of the sample's assigned reads |
| `results/consolidated_output/otu_table_counts.tsv` | Species × samples: estimated read counts, including the `Unassigned` reads |
| `results/consolidated_output/tax_table.tsv`, `sample_metadata.tsv` | The lineage of each species, and per sample the species detected and the assigned and unassigned reads |
| `results/sample_outputs/<sample>_rel-abundance.tsv` | One table per sample, with each species' lineage, counts and the warning columns explained below |

The `Below detection` and `Unassigned` rows are 0 in `otu_table_abundance.tsv` by design: abundances are shares of the reads that were assigned to a reported species. How many reads were not assigned is in `otu_table_counts.tsv` and `sample_metadata.tsv`. [What the numbers mean](#what-the-numbers-mean) explains the columns that flag uncertain species, and [Output](#output) lists every file.

## How it works

| Step | What it does | Tool | Origin |
|------|--------------|------|--------|
| Filter | Keeps reads of 500-2,000 bases (`--min_length`, `--max_length`); no quality filter unless `--min_quality` is set | [fastplong](https://github.com/OpenGene/fastplong) 0.7.0 | external tool |
| Index and align | Builds a minimap2 index of the database once, then aligns each read to it (`-ax map-ont -N 3 -p 0.9`): its best match and up to 3 other candidate sequences | [minimap2](https://github.com/lh3/minimap2) 2.31, samtools 1.24 | external tools |
| Score | Scores each read against each candidate species from its alignment (the CIGAR string), with error rates learned from the sample's own alignments | NanoVI (`bin/likelihood.py`) | adapted from [Emu](https://github.com/treangenlab/emu)'s read likelihood |
| Estimate | Estimates each species' share of the reads with a Bayesian model (a Dirichlet prior, fitted by variational inference), drops species below the detection floor and refits | NanoVI (`bin/variational_inference_bayesian.py`) | NanoVI's own; `--estimator em` runs Emu's EM loop instead |
| Write | Writes the per-sample and whole-run tables, and the optional outputs | NanoVI (`bin/`) | NanoVI's own; the table layout follows Emu |
| Report | With `--report`: an HTML site and PDF with figures | NanoVI (`bin/nanovi_report.py`) | NanoVI's own |

Reads are filtered by length only: the 2026 benchmark found that a mean-quality filter changed little on real data and distorted the composition where it removed reads ([methods_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/methods_2026.md) §7, in [NanoVI-paper](https://github.com/microbialds/NanoVI-paper)). The model and the evidence behind each default are in [docs/methods.md](docs/methods.md); every step of the pipeline, the parameter that switches it on and its container image are in [docs/workflow.md](docs/workflow.md).

NanoVI has four commands, chosen with `--cmd`:

| Command | What it does |
|---------|--------------|
| `abundance` (the default) | Estimate species abundances from reads |
| `build-database` | Build a reference database ([Reference databases](#reference-databases)) |
| `collapse-taxonomy` | Sum a per-sample table to a higher rank, e.g. genus ([Other commands](#other-commands)) |
| `combine-outputs` | Merge per-sample tables into one table at a chosen rank ([Other commands](#other-commands)) |

## Input: samplesheet or folder

`--input` takes either a samplesheet or a folder.

- **Use a samplesheet** (recommended) when you want to choose the sample names, when the files are in different folders or on remote storage, or when their file names are not the names you want in the results.
- **Use a folder** when all the FASTQ files are in one folder (subfolders are searched too) and each file name, without its extension, is the sample name you want.

### Samplesheet

A CSV file with `sample` and `fastq` columns:

```csv
sample,fastq
sample1,/path/to/sample1.fastq.gz
sample2,reads/sample2.fastq.gz
sample3,s3://my-bucket/run1/sample3.fastq.gz
```

It is checked against `assets/schema_input.json` before any step runs: sample IDs must be unique, start with a letter or digit and may use letters, digits, `.`, `_` and `-` (as in folder mode); each `fastq` must exist and end in `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz`, with no spaces in the path. Paths may be absolute, relative to the directory you launch Nextflow from (not to the samplesheet), or remote (`s3://`, `gs://`, `az://`, `https://`). Fields may be quoted. Every line after the header is a sample: lines starting with `#` are not skipped, and a sheet with no rows stops the run.

### Folder

A path to a folder containing FASTQ files (`.fastq`, `.fq`, optionally gzipped), searched recursively. The sample ID is the file name without that extension, so `S1.run1.fastq.gz` becomes `S1.run1`. If two files give the same ID (for example `S1.fastq` and `S1.fq`, or the same name in two subfolders), the pipeline stops and lists them; rename the files or use a samplesheet. IDs must start with a letter or digit and may use letters, digits, `.`, `_` and `-`.

## Parameters

Parameters start with two dashes (`--input`); Nextflow's own options, such as `-profile`, `-r` and `-resume`, with one. Parameters are checked against `nextflow_schema.json` at the start of every run (types, allowed values, ranges, and that input paths exist), and values given on the command line are converted to the declared types, so `--keep_counts false` turns the counts off. A parameter the pipeline does not know stops the run: a misspelling such as `--output-dir` gives "Parameter `outputDir` was specified on the command line or params file but is not declared". In the tables, the second column is the value used when you do not set the parameter.

### Required

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--input` | none: you must set it | Samplesheet CSV or folder of FASTQ files ([Input](#input-samplesheet-or-folder)) |
| `--db` | none: you must set it | Database directory containing `species_taxid.fasta` and, unless `--taxonomy_tsv` is given, `taxonomy.tsv` ([Reference databases](#reference-databases)). Never a `.mmi` file: pass a prebuilt index with `--db_index` |

### Common options

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--cmd` | `abundance` | The command to run: `abundance`, `build-database`, `collapse-taxonomy` or `combine-outputs` ([How it works](#how-it-works)) |
| `--output_dir` | `results/` | Folder for the results. A second run into the same folder replaces the files it writes; files that only an option writes stay from the first run when the second does not use that option ([Troubleshooting](#troubleshooting)) |
| `--keep_counts` | `true` | Add the `estimated_counts` column to the per-sample tables and write `otu_table_counts.tsv` |
| `--binomial_names` | `false` | Cut species names in the output tables to genus and epithet (`Escherichia coli O157` → `Escherichia coli`; `Candidatus` names keep three words) and sum the species that then share a name (their `shared_read_fraction` counts only the reads split with species of other rows). Off: names are kept as the database's `taxonomy.tsv` gives them |
| `--credible_intervals` | `false` | Add `ci_lower` and `ci_upper`, 95% sampling intervals of each species' abundance (`vb` only; [what the numbers mean](#what-the-numbers-mean)) |
| `--ambiguity_groups` | `false` | Also write `sample_outputs/<sample>_groups.tsv`: the rows of the abundance table grouped where the sample cannot separate them, with their summed abundance ([Optional outputs](#optional-outputs)). The abundance table does not change |
| `--cell_abundance` | `false` | Also write `sample_outputs/<sample>_cell-abundance.tsv`: each species' abundance divided by its 16S copy number and renormalised, an estimate of its share of cells ([Optional outputs](#optional-outputs)). Needs a database built with copy numbers; the abundance table does not change |
| `--read_assignments` | `false` | Also write `sample_outputs/<sample>_read_assignments.tsv.gz`: for every read that passed the length filter, the rows of the abundance table that take at least 0.01 of it and their probability, or `Below detection` or `Unassigned` ([Optional outputs](#optional-outputs)) |
| `--report` | `false` | Also write `report/`: an HTML site that works offline (`report/index.html`), `nanovi_report.pdf` and the figures as SVG, PDF, 600 dpi PNG and CSV ([Optional outputs](#optional-outputs)). The other outputs do not change |
| `--report_metadata` | not used | Sample metadata for the report, TSV or CSV with a `sample_id` column and optionally `group`: samples are drawn by group (in the order the file first names each group), then in natural order (S2 before S10); other columns are listed in the report |
| `--index_cache_dir` | `<output_dir>/index_cache` | Where the alignment index of the database is kept. NanoVI creates it; set it only to reuse one index across projects (point every project at the same folder) |
| `--keep_files` | `false` | Also publish the alignments (`bam/<sample>_vi_alignments.bam`) and the filtered reads (`filtered/<sample>_filtered.fastq.gz`) |

### Advanced: reads and alignment

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--min_length` | `500` | Shortest read kept, in bases |
| `--max_length` | `2000` | Longest read kept, in bases |
| `--min_quality` | not set: no quality filter | Minimum mean read quality (Phred). Set: fastplong's quality filter drops reads with a lower mean quality (its per-base rule is switched off; its limit of 5 N bases per read still applies) |
| `--taxonomy_tsv` | `<db>/taxonomy.tsv` | Taxonomy table, if not the database's own |
| `--db_index` | not set: NanoVI builds the index | Prebuilt minimap2 index (`.mmi`) of the database's `species_taxid.fasta`; the index step is skipped. Build it with the pipeline's minimap2 (2.31) and the same `-x <type> -k <kmer_size>`, or take it from the index cache of an earlier run |
| `--kmer_size` | `21` | K-mer size for the minimap2 index and mapping (`-k`, 1-28) |
| `--N` | `3` | Most secondary alignments per read |
| `--K` | `500M` | Minibatch size for minimap2 mapping (minimap2 `-K`): the number of read **bases** (not bytes) loaded and mapped per batch, a number optionally with `K`, `M` or `G`. 500M (minimap2's default) is 500 million bases, about 330,000 reads of 1.5 kb; larger batches use more memory |
| `--type` | `map-ont` | Minimap2 preset (`map-ont`, `map-pb`, `sr`), used for the index and the alignment |
| `--split_prefix` | `temp` | Prefix of minimap2's temporary files when it splits a large index (`--split-prefix`) |

**The index cache.** The index file's name records the database (the first 12 digits of the FASTA's MD5), `--type`, `--kmer_size` and the minimap2 version, e.g. `species_taxid.3f9c0a12b7de.map-ont.k21.minimap2-2.31.mmi`, so a later run with the same settings reuses the index instead of rebuilding it. The MD5 is read once per FASTA and kept in `md5_cache.tsv` in the same directory, keyed by the FASTA's path, size and modification time (and, for a local file, its inode and status-change time, so a replaced file is read again), so later launches do not read the database again (delete the file to make the next run recompute it). Point several runs at one directory to share it, but let the first run that builds an index finish before starting others: a run reuses any file with the right name, and nothing checks that it is complete. Index-time settings (k-mer size, window, homopolymer compression) are fixed in an index, and an index from another minimap2 version may not load.

### Advanced: abundance estimation

| Parameter | Default (used when you don't set it) | Estimator | Description |
|-----------|---------|-----------|-------------|
| `--estimator` | `vb` | | `vb`: Bayesian variational inference (Dirichlet prior, CAVI updates). `em`: the expectation-maximization loop of Emu, kept as a baseline |
| `--min_abundance` | `0.001` | both | Lowest abundance reported (see the detection floor below) |
| `--vi_tolerance` | `1e-7` | vb | CAVI stops when no species' posterior mean abundance changes by more than this between iterations |
| `--vi_max_iterations` | `1000` | vb | Most CAVI iterations in each pruning round |
| `--alpha_prior` | `1.0` | vb | Total concentration of the symmetric Dirichlet prior; each of the S candidate species gets `alpha_prior / S` |
| `--vi_elbo_threshold` | `1e-4` | em | The EM loop stops when the log-likelihood rises by less than this. Its inner loop keeps Emu's settings (at most 20 iterations, total change in frequencies below 5e-4) |

### Advanced: computing resources

The Python steps ask for one CPU each and run on one core (`nextflow.config` sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` to 1). minimap2 asks for 8 CPUs; when a task is killed for memory or time it is retried (twice) with more memory and time but the same CPUs. No request goes above these caps (`process.resourceLimits`):

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--max_cpus` | `16` | Most CPUs any task asks for |
| `--max_memory` | `128.GB` | Most memory any task asks for |
| `--max_time` | `240.h` | Most time any task asks for |

Give the caps on the command line or in a `-params-file`. In a custom `-c` config, `params.max_*` comes too late to change them: set `process.resourceLimits = [cpus: 32, memory: 256.GB, time: 48.h]` there instead. The same holds for `--keep_files`: a `params.keep_files` in a `-c` config is ignored. On a cluster, also review the per-step requests in `conf/base.config` against your partitions.

The parameters of `--cmd build-database` are under [Reference databases](#reference-databases), and those of `collapse-taxonomy` and `combine-outputs` under [Other commands](#other-commands).

## What the numbers mean

**Abundances are shares of 16S reads.** Each species' `abundance` is its share of the reads assigned to reported species; the shares add up to 1 in every sample. A species with more copies of the 16S gene per genome gives more reads per cell, so these are not shares of cells: `--cell_abundance` estimates those ([Optional outputs](#optional-outputs)). The `Below detection` and `Unassigned` rows have abundance 0 and carry reads only in `estimated_counts`.

**Detection floor.** R is the number of reads with at least one scored alignment. With `vb`, species whose posterior mean abundance falls below the pruning threshold are removed and the rest refitted, until none is below it: the threshold is max(1/R, `min_abundance`) up to 1,000 reads and max(10/R, `min_abundance`) above. With the default, the floor is 1/R up to 1,000 reads, 10/R up to 10,000 reads and 0.1% above that: nothing below 0.1% is reported at any depth (Emu's floor of 10/R would be 0.01% at 100,000 reads). With `em`, species are pruned below 1/R (10/R above 1,000 reads), as in Emu, and if `min_abundance` is higher the species at or above it are refitted and reported. With either estimator, reads whose candidate species were all pruned are reported in the `Below detection` row; reads that align to nothing (or not well enough to be scored) are `Unassigned`.

**Shared reads and intervals.** Every per-sample table has a `shared_read_fraction` column: the share of the species' expected reads whose responsibility is split with another reported species, that is, where the other species get at least 1% of the read. It is counted per row of the table: in a row that sums several tax IDs (with `--binomial_names`, or a name repeated in a custom taxonomy), reads split only among those tax IDs are not shared. A high value means the species' abundance rests on reads that another species explains as well; in the 2026 benchmark it flagged the estimates that missed the truth (AUC 0.75-0.84; [methods_2026.md](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/methods_2026.md) §8). With `--credible_intervals` (`vb` only) the table also gets `ci_lower` and `ci_upper`, the 2.5% and 97.5% quantiles of the species' mean-field posterior marginal, Beta(α_s, Σα − α_s) over the reported species (tax IDs that share a species name add their α). These are sampling intervals, not calibrated uncertainty: they shrink as 1/√reads while the error from assigning reads between close relatives does not, so they covered the true abundance in 93%, 81% and 38% of cases at 1,000, 10,000 and 50,000 reads (42% for species with a close relative at 10,000), and the exact posterior (collapsed Gibbs) and a read bootstrap gave the same intervals. The columns are left empty on the `Below detection` and `Unassigned` rows, and `collapse-taxonomy`, `combine-outputs` and the consolidated tables leave them out.

**Species with identical 16S.** The per-sample table's `identical_16s_species` column names, for each reported species, the species of the same genus in the database that share an identical 16S sequence with it (`; `-separated, empty for most species). 16S cannot tell such species apart: reads that match both score the same, so the reported split between them is not evidence of either. The list comes from the database's `identical_16s.tsv`, which `--cmd build-database` writes; for a database without it, the pipeline builds it once from `species_taxid.fasta` and `taxonomy.tsv` and keeps it beside the cached indexes in `--index_cache_dir`, under a name made of the MD5s of the FASTA, the taxonomy and `bin/identical_16s.py`. An `identical_16s.tsv` in `--db` is used as it is when the taxonomy is `<db>/taxonomy.tsv` (with another `--taxonomy_tsv` the table is built from the FASTA and that taxonomy); delete it to have the pipeline build the table again. With `--db_index` and no FASTA in `--db` the column is left empty.

## Optional outputs

**Per-read assignments.** With `--read_assignments` each sample also gets `<sample>_read_assignments.tsv.gz`, a long table of `read_id`, `species` and `probability`. It lists every read that passed the length filter, with one line for each row of the abundance table that takes at least 0.01 of the read. The probability is the read's responsibility in the final fit (φ for `vb`, the EM posterior for `em`), summed over the row's tax IDs, so `--binomial_names` merges rows here as in the table. A read's probability on species that are not reported is its `Below detection` line, and reads with no scored alignment are `Unassigned` with probability 1. Use it to see which reads a species' abundance rests on, or which reads are split between it and another species (`shared_read_fraction` uses the same 0.01). The file is not written by default, and only then does the scoring step keep the read names ([docs/output.md](docs/output.md)).

**Ambiguity groups.** Species whose 16S genes are identical or nearly so give reads that fit both, and the estimator splits such reads by abundance, so the abundance table can report, say, *E. coli* at 0.26 and *S. flexneri* at 0.04 where the truth is 0.25 and 0.05, or one species of a pair at 0. With `--ambiguity_groups` each sample also gets `<sample>_groups.tsv`, which reports such rows together, with the sum the reads support. Two rows are linked when they are in the same family and at least 0.25 of either one's expected reads sit on reads that the other row takes at least 0.01 of, and a group is a set of linked rows; every row is in one group, alone if it links to nothing (as is a row with no family in `taxonomy.tsv`). Each group also lists `identical_16s_relatives`: the species that share an identical 16S with a member in the database and are not reported, which the group may also contain. Groups come from each sample's reads, so they can differ between samples; the consolidated tables, `collapse-taxonomy` and `combine-outputs` use the species tables ([docs/output.md](docs/output.md)).

**Cell-level abundance.** The abundance table gives each species' share of 16S reads, which follows its share of 16S genes, not of cells: a species with ten copies of the gene per genome gives about ten times the reads per cell of a species with one. With `--cell_abundance` each sample also gets `<sample>_cell-abundance.tsv`, where each species' abundance is divided by its copy number from the database's `copy_numbers.tsv` and the result renormalised. A copy number is the median count over the species' complete genomes in GTDB; a species without one takes the median of its genus, or else of its family, and the `copy_number_level` column says which (`species`, `genus` or `family`). Species with none of the three, and the `Below detection` and `Unassigned` reads, have no cell-level share and are left out; the file's first two lines (starting with `#`) say how many reads and how much of the abundance that is. A row that sums several tax IDs (`--binomial_names`), and with `--ambiguity_groups` a group, takes the cell-weighted mean of its members' copy numbers, so its cell share is the sum of theirs. A database built without copy numbers stops the run with an error, as does a `--taxonomy_tsv` other than the database's own (the copy numbers are keyed by its tax IDs); build it with `--gtdb_metadata`, or with `--copy_numbers` / `--rrndb_pantaxa` ([Reference databases](#reference-databases)). The consolidated tables, `collapse-taxonomy` and `combine-outputs` use the species tables ([docs/output.md](docs/output.md)). In the 2026 benchmark ([NanoVI-paper `copy_number_2026.md`](https://github.com/microbialds/NanoVI-paper/blob/main/docs/benchmarks/copy_number_2026.md)), on simulated communities whose reads follow cells × copies, the cell table's species L1 against the cell composition was 0.130 where the read-level table's was 0.343; on 12 Zymo mocks 0.367 against 0.438. A copy number filled from the genus median instead of the species' own value added 0.106 and 0.033. The correction cannot fix species the read-level estimate misses or splits wrongly.

**Run report.** With `--report` a last step reads what the run published (the species tables, the consolidated tables, fastplong's reports, `software_versions.yml`, and the group, cell-level and per-read tables when their options ran) and writes `report/`. Its sections are the read fate of each sample (raw reads, reads that passed the length filter, and how many of those went to reported species, `Below detection` and `Unassigned`), the composition at species and genus (the union of each sample's top 5, at most 20, the rest as Other; a heatmap of the top 30 species above 14 samples), the flagged species (a `shared_read_fraction` of 0.5 or more in a sample, or a non-empty `identical_16s_species`; the cutoff is a display rule, not a test), the ambiguity groups and the cell-level table when present, a summary of the per-read file, and the methods: the run's parameters, every tool and version with its reference (Emu credited) and how to read the tables. The HTML pages load nothing from the network: the fonts and the Vega libraries of the interactive read-fate and composition figures are in `report/assets/`. The figures follow journal widths (89 or 183 mm) with Inter embedded; a figure of more than 60 samples is split into parts. The report runs in its own image, `ghcr.io/microbialds/nanovi-report:<version>` (Typst, Jinja, matplotlib, Altair; `containers/report/`), so the Python image does not change. If Typst fails, the run warns and writes the site without the PDF. The report carries the day the run started, and two renders of the same results give the same files. `bin/nanovi_report.py --results <output_dir>` renders one from a finished run's directory in that image. [docs/output.md](docs/output.md) lists its files.

## Output

`--cmd abundance` writes into `--output_dir`:

| Path | Description |
|------|-------------|
| `sample_outputs/<sample>_rel-abundance.tsv` | One table per sample: a row per species (abundance, lineage, `estimated_counts`, `identical_16s_species`, `shared_read_fraction`, and `ci_lower`/`ci_upper` with `--credible_intervals`), then the `Below detection` and `Unassigned` rows |
| `sample_outputs/<sample>_read_assignments.tsv.gz` | Only with `--read_assignments`: every read's table rows and probabilities (above) |
| `sample_outputs/<sample>_groups.tsv` | Only with `--ambiguity_groups`: the species rows grouped where the sample cannot separate them, with summed abundance and counts (above); with `--cell_abundance` also their cell-level share and copy number |
| `sample_outputs/<sample>_cell-abundance.tsv` | Only with `--cell_abundance`: each species' abundance divided by its 16S copy number and renormalised, with the copy number and the level it was filled from (above) |
| `consolidated_output/otu_table_abundance.tsv` | Species × samples abundance matrix (phyloseq `otu_table`); `Below detection` and `Unassigned` are 0 here (above) |
| `consolidated_output/otu_table_counts.tsv` | Species × samples `estimated_counts` matrix, with the `Below detection` and `Unassigned` reads; not written when no sample has counts (`--keep_counts false`) |
| `consolidated_output/tax_table.tsv` | Species × ranks (genus to superkingdom; phyloseq `tax_table`) |
| `consolidated_output/sample_metadata.tsv` | Per sample: total abundance, species detected, assigned and unassigned reads |
| `qc/fastplong/<sample>_report.html.gz`, `qc/fastplong/<sample>_report.json.gz` | fastplong's report: reads and bases before and after the length filter |
| `filtered/<sample>_filtered.fastq.gz` | Reads that passed the length filter (only with `--keep_files`) |
| `bam/<sample>_vi_alignments.bam` | minimap2 alignments, primary and up to `--N` secondary per read (only with `--keep_files`) |
| `index_cache/` | The minimap2 index and, for a database without `identical_16s.tsv`, the identical-16S table the pipeline builds, each with a `.versions.yml`; reused by later runs (elsewhere with `--index_cache_dir`; no index with `--db_index`) |
| `pipeline_info/timeline.html`, `report.html`, `trace.txt` | Nextflow's timeline, report and trace (with each task's requested CPUs, memory and time), replaced by the next run into the same folder |
| `pipeline_info/software_versions.yml` | Versions of fastplong, minimap2, samtools, Python and its packages, NanoVI and Nextflow used in the run |
| `report/` | Only with `--report`: the run report, `index.html` (an offline site), `nanovi_report.pdf`, `figures/` and the report's own `versions.yml` (above) |

Every column, the consolidated tables and the database files are described in [docs/output.md](docs/output.md); the model behind the numbers, with the benchmark evidence for each default and option, is in [docs/methods.md](docs/methods.md), and which step writes each file in [docs/workflow.md](docs/workflow.md).

## Reference databases

A NanoVI database is a directory with `species_taxid.fasta` (the 16S sequences, one tax ID per species) and `taxonomy.tsv` (each species' lineage), plus `identical_16s.tsv`, `cross_genus_identical.tsv` and, for `--cell_abundance`, `copy_numbers.tsv`. You pass the directory as `--db`.

**GTDB is the recommended database.** NanoVI's defaults were chosen and benchmarked on GTDB release r226 ([docs/methods.md](docs/methods.md)), built as in [step 3](#3-build-the-reference-database-once). GTDB's `latest` folder holds the newest release (r232 since April 2026) with the release number left out of the file names (`ssu_all.fna.gz`, `bac120_metadata.tsv.gz`); each release's own folder, such as `releases/release232/232.0/`, has the numbered names (`ssu_all_r232.fna.gz`). r232 also builds and runs with NanoVI 1.2.0, but has not been benchmarked, and it changes results: it adds species such as *ECMA0423 sp047199055*, a new Enterobacteriaceae genus whose 5,575 database sequences include 634 identical to one of *Escherichia coli*, so on NanoVI's test reads r232 gives *E. coli* 0.17 and *ECMA0423 sp047199055* 0.11 where r226 gives *E. coli* 0.26. Use r226 for results comparable with NanoVI's benchmarks and publications. Other databases can be built from a FASTA and a taxonomy (below).

### From a GTDB SSU file

`--cmd build-database --gtdb_ssu <file>` runs `bin/build_gtdb_db.py` in the pipeline's container. It reads each sequence's lineage from GTDB's FASTA headers, assigns **one tax ID per species**, keeps identical sequences of a species once, and writes `<output_dir>/<db_name>/`.

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--gtdb_ssu` | none: you must set it | GTDB's SSU FASTA of bacteria and archaea (`ssu_all_r226.fna.gz`) |
| `--db_name` | `db_custom` | Name of the database directory, also written into the FASTA headers (e.g. `db_gtdb_r226`) |
| `--gtdb_metadata` | not used: no copy numbers | The release's `bac120_metadata` and `ar53_metadata` files, comma-separated or as a glob (`gtdb_r226/{bac120,ar53}_metadata_r226.tsv.gz`): write `copy_numbers.tsv`, each species' 16S copy number as the median `ssu_count` (16S genes found in the genome) over its `Complete Genome` assemblies; draft assemblies are not used, because their 16S copies collapse into one or a few contigs (in r226, 72% of drafts count 0 or 1). A species without a complete genome takes the median of its genus, or else of its family. With `--merge_gtdb_suffixes` the merged species pool their genomes |
| `--copy_numbers` | not used | Instead of `--gtdb_metadata`, a TSV with `species` and `copy_number` columns, filled the same way |
| `--gtdb_min_length` | `900` | Shortest 16S sequence taken into the database, in bp (not `--min_length`, which filters reads) |
| `--merge_gtdb_suffixes` | `false` | Strip GTDB's clade suffixes (`Listeria monocytogenes_B` → `Listeria monocytogenes`, `Pseudomonas_E` → `Pseudomonas`), merging the species and genera GTDB separates with them into one tax ID, and write `merged_species.tsv`. Off: names are kept as GTDB writes them |
| `--cross_genus_rule` | `keep` | What to do with the sequences found in species of more than one genus (`cross_genus_identical.tsv`): `keep` them; `remove` every copy; keep each only in the genus with the most species carrying it (`best-genus`); or only in the genus whose entries stand for the most input sequences, about one per GTDB genome (`most-copies`; a tie goes to the genus with more entries in the database). A tie that the rule cannot break removes every copy. Species left with no sequence stay in `taxonomy.tsv`; `identical_16s.tsv` describes the cleaned database. Writes `cross_genus_removed.tsv` and `cross_genus_emptied.tsv`. Off by default: on the benchmark the rules remove a few false-positive species but make the real mocks' abundances slightly worse ([docs/output.md](docs/output.md)) |

The same build runs without Nextflow from a clone of this repository, with the Python packages in `containers/requirements.txt`: `python3 bin/build_gtdb_db.py --ssu ssu_all_r226.fna.gz --db-name db_gtdb_r226 --output-dir db_gtdb_r226 --gtdb-metadata bac120_metadata_r226.tsv.gz ar53_metadata_r226.tsv.gz`. The script spells its options with dashes and also accepts the pipeline's underscores (`--output_dir`).

**Output files:**

| File | Description |
|------|-------------|
| `species_taxid.fasta` | Reference FASTA — one entry per unique 16S sequence, header: `taxid:db_name:n` |
| `taxonomy.tsv` | Taxonomy table — one row per species with full lineage |
| `identical_16s.tsv` | Species that share an identical 16S sequence with another species of their genus (genera compared without GTDB's suffixes), one row each: its group (the species of the genus linked through shared sequences) and the species it shares a sequence with. 16S cannot tell these species apart; the pipeline names them in the `identical_16s_species` column |
| `cross_genus_identical.tsv` | Sequences found in species of more than one genus, with their species, genera, families and phyla. In GTDB most are contaminant contigs in an assembly, but real strains carry some of them (1,360 sequences in r226). With the default `--cross_genus_rule keep` they stay in the database, and `identical_16s.tsv` still groups the species of one genus that carry such a sequence; only the link across genera is left out |
| `cross_genus_removed.tsv` | With a `--cross_genus_rule` other than `keep`: one row per database entry the rule removed (sequence ID, tax ID, species, genus, the sequence's MD5, the rule and the genus the sequence was kept in, empty when none was) |
| `cross_genus_emptied.tsv` | With a `--cross_genus_rule` other than `keep`: the species the rule left with no sequence (tax ID, species, genus, entries removed). They stay in `taxonomy.tsv` but can no longer be reported |
| `copy_numbers.tsv` | With `--gtdb_metadata` or `--copy_numbers`: one row per species of `taxonomy.tsv` with its 16S copy number, the level it comes from (`species`, `genus` or `family` median, or `none`), the number of complete genomes behind a species value and the source. GTDB r226: 14% of the default build's 84,596 species have their own value, 37% take a genus median, 26% a family median and 24% none. `--cell_abundance` reads it |
| `merged_species.tsv` | With `--merge_gtdb_suffixes` only: each name that merges two or more GTDB species with a sequence in the database, with those species and their genera, the ranks at which their lineages differ, and the GTDB species whose lineage `taxonomy.tsv` keeps (the first in the input) |

### From your own sequences and taxonomy

For a database other than GTDB, `--cmd build-database` takes a FASTA, a table linking each sequence to a tax ID, and either a taxonomy table or an NCBI taxonomy dump:

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--db_name` | `db_custom` | Name for the custom database |
| `--sequences` | none: you must set it | Input FASTA with reference sequences |
| `--seq2tax` | none: you must set it | Sequence-to-taxonomy mapping TSV |
| `--taxonomy_list` | set this or `--ncbi_taxonomy` | Taxonomy TSV with a `tax_id` column and one column per rank |
| `--ncbi_taxonomy` | set this or `--taxonomy_list` | NCBI taxonomy dump directory (`names.dmp`, `nodes.dmp`) |
| `--cross_genus_rule` | `keep` | As above |
| `--copy_numbers` | not used | TSV with `species` and `copy_number` columns (16S copies per genome): write `copy_numbers.tsv` for `--cell_abundance`, matched on the species names of the built taxonomy; species without a value take the median of their genus, or else family. At most one of `--copy_numbers` and `--rrndb_pantaxa` |
| `--rrndb_pantaxa` | not used | [rrnDB](https://rrndb.umms.med.umich.edu/)'s pan-taxa statistics for NCBI taxonomy (`rrnDB-<version>_pantaxa_stats_NCBI.tsv`): write `copy_numbers.tsv` from the median of each species' NCBI taxid, filled the same way. Needs `--ncbi_taxonomy` |

It writes `<output_dir>/<db_name>/` with `species_taxid.fasta`, `taxonomy.tsv`, `identical_16s.tsv` and `cross_genus_identical.tsv`, plus `cross_genus_removed.tsv` and `cross_genus_emptied.tsv` with a `--cross_genus_rule` other than `keep`, and `copy_numbers.tsv` with `--copy_numbers` or `--rrndb_pantaxa` (see [docs/output.md](docs/output.md)); pass that directory as `--db`. NanoVI's defaults were not benchmarked on other databases.

## Other commands

| Parameter | Default (used when you don't set it) | Description |
|-----------|---------|-------------|
| `--input_tsv` | none: `collapse-taxonomy` needs it | Abundance table for taxonomy collapse |
| `--input_dir` | none: `combine-outputs` needs it | Directory with per-sample tables for combining |
| `--rank` | none: both need it | Taxonomic rank: `species`, `genus`, `family`, `order`, `class`, `phylum`, `superkingdom` |

`--cmd collapse-taxonomy` writes `<output_dir>/<input name>-<rank>.tsv` (for `S1_rel-abundance.tsv`, `S1_rel-abundance-genus.tsv`), and `--cmd combine-outputs` writes `<output_dir>/vi-combined-<rank>.tsv`, one column per `*_rel-abundance.tsv` in `--input_dir`; a directory without such tables stops the command with an error.

## Profiles: Docker, Singularity, conda and clusters

A profile tells Nextflow how to run the tools. For a first run, `-profile docker` is all you need.

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (also used if no profile is given) |
| `singularity` | Run with Singularity/Apptainer containers, for clusters without Docker |
| `conda` | Run with Conda environments: fastplong 0.7.0 and minimap2 2.31 (with samtools 1.24) from bioconda, and `containers/environment.yml` (the image's Python and package versions) for the Python steps; with `--report`, `containers/report/environment.yml`. Creating the environments with conda's classic solver can take a long time; a recent conda (libmamba solver) or `-c` with `conda.useMicromamba = true` is much faster |
| `slurm` | Submit jobs to the SLURM scheduler (uses Singularity) |
| `sge` | Submit jobs to the SGE scheduler (uses Singularity) |
| `test` | The bundled test data (`tests/test_data/e2e`) with reduced resources; it sets no container engine, so combine it with one: `-profile test,docker` |

Before running on a cluster, review the resource requests (CPUs, memory, time) in `conf/base.config` against your partitions; `--max_cpus`, `--max_memory` and `--max_time` cap every request, retries included ([Parameters](#advanced-computing-resources)).

## Troubleshooting

- **"Parameter `outputDir` was specified on the command line or params file but is not declared".** The pipeline's parameters are written with underscores: `--output_dir`, not `--output-dir` (the same for any other parameter). Only `bin/build_gtdb_db.py`, when run on its own, uses dashes, and it accepts underscores too.
- **Running again into the same `--output_dir`.** The new run replaces the results and the files in `pipeline_info/`. Add `-resume` to reuse the steps that already finished. Nothing is deleted: files that only an option writes (`<sample>_groups.tsv`, `<sample>_cell-abundance.tsv`, `<sample>_read_assignments.tsv.gz`, `report/`, `bam/`, `filtered/`) stay from the earlier run when the new run does not use that option, and so do the tables of samples no longer in the input; a database rebuilt into the same `<output_dir>/<db_name>/` keeps the earlier build's files that the new build does not write (`copy_numbers.tsv`, `merged_species.tsv`, `cross_genus_removed.tsv`, `cross_genus_emptied.tsv`). When you change the options or the samples, use a new `--output_dir` (or delete the old one). Before NanoVI 1.2.0, a second run warned "Failed to create trace file … already exists" and "Failed to render execution report"; the results were written, only those Nextflow files were not.
- **"Cannot connect to the Docker daemon" or "permission denied" with Docker.** Docker is not running, or your user may not use it: see [step 1](#1-install-the-requirements). On a cluster, use `-profile singularity`.
- **The run report's image cannot be pulled.** `--report` uses `ghcr.io/microbialds/nanovi-report`; if it cannot be pulled, run without `--report` or build it ([docs/development.md](docs/development.md)).
- **A sample gets only `Unassigned`.** None of its reads aligned to the database (for example, they are not 16S, or are shorter than `--min_length`); the table says so instead of failing the run.
- **Where is my FASTQ?** Paths in a samplesheet are relative to the folder you launch Nextflow from, not to the samplesheet; absolute paths avoid the question.

Developers: building the containers, the tests and the release steps are in [docs/development.md](docs/development.md).

## Credits

NanoVI builds on Emu (Curry KD, Wang Q, Nute MG, et al. Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data. *Nature Methods* 19, 845–853 (2022). [doi:10.1038/s41592-022-01520-4](https://doi.org/10.1038/s41592-022-01520-4); [github.com/treangenlab/emu](https://github.com/treangenlab/emu)). Parts of NanoVI's Python code are adapted from Emu's source:

- **The CIGAR-based read score.** The log-likelihood of a read under each candidate species is Emu's: insertion, deletion, soft-clip and mismatch rates are estimated from the sample's own alignments, and each alignment is scored from its CIGAR counts, scaled by the length of the read's longest alignment (`bin/likelihood.py`).
- **The EM loop.** The expectation-maximization loop with pruning of low-abundance species in `bin/em_estimator.py` (`--estimator em`) follows Emu's. The default estimator, the Bayesian variational inference in `bin/variational_inference_bayesian.py`, is NanoVI's own: it places a Dirichlet prior on the species frequencies and estimates their posterior instead of a point estimate.
- **The database layout.** A `species_taxid.fasta` with `taxid:db_name:n` headers plus a `taxonomy.tsv` lineage table, the build-database helpers, the per-sample output table and the collapse-taxonomy and combine-outputs commands follow Emu (`bin/taxonomy.py`, `bin/main.py`).

Emu is distributed under the MIT License, Copyright (c) 2021 Kristen Curry, Qi Wang, Michael G. Nute, Alona Tyshaieva, Qinglong Wu, Tor Savidge, Alexander Dilthey, Todd J. Treangen. Its full license text is in [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES), and each file with adapted code names its origin in its header. If you use NanoVI, please cite Emu as well.

## Citation

If you use NanoVI in your research, please cite the preprint:

> Curiqueo C, Fuentes-Santander F, Ugalde JA. NanoVI: a Bayesian variational inference Nextflow pipeline for species-level taxonomic classification from full-length 16S rRNA Nanopore reads. *bioRxiv* 2026. [doi:10.64898/2026.03.07.710315](https://doi.org/10.64898/2026.03.07.710315)

and Emu, whose read likelihood NanoVI uses (see Credits). The citation metadata are also in [CITATION.cff](CITATION.cff), which GitHub shows as "Cite this repository". This section will give the journal article once it is published.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details. Code adapted from Emu remains under Emu's MIT License; see [THIRD_PARTY_LICENSES](THIRD_PARTY_LICENSES).
