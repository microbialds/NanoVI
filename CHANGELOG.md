# Changelog

All notable changes to the NanoVI pipeline will be documented in this file.

## [Unreleased]

### Added

- Clinical rerun for the manuscript (`benchmarks/clinical.py`, `docs/benchmarks/clinical_2026.md`, `docs/benchmarks/data/clinical_*.tsv`): NanoVI on GTDB r226 and on Emu's database, Emu 3.6.2 as shipped and Emu on NanoVI's alignments, on three ONT 16S cohorts (149 vaginal and stool samples). On the same alignments NanoVI and Emu give the same profiles (median species Bray-Curtis 0.004 or less); profiles on GTDB r226 differ from those on Emu's database through what each database holds. The manuscript figures add Figure 8 (W10, D34).
- `CITATION.cff` with the preprint (bioRxiv, doi:10.64898/2026.03.07.710315) as the preferred citation and Emu as a reference; the README's Citation section cites the preprint instead of a placeholder (W11).
- Every process writes an nf-core style `versions.yml`, and each run collects them into `pipeline_info/software_versions.yml`: the fastplong and minimap2 versions, the Python and package versions of every Python step, and the NanoVI and Nextflow versions (F33).
- Unit tests for the code the pipeline runs: the Bayesian VI (every update against a hand or Monte Carlo calculation, the ELBO, the stopping and pruning rules, known proportions recovered from synthetic reads, and agreement within 1e-9 with the earlier dense implementation), the likelihood on small SAM files, the output and consolidated tables, the database reports and every `bin/` script; 227 tests, run in CI on Python 3.11 and 3.13 (F36).
- `docs/methods.md`: the model in equations (Emu's read score, the Dirichlet-categorical model, the mean-field updates, pruning, `shared_read_fraction` and the sampling intervals), what differs from Emu and why, each default with its benchmark evidence, the database construction and the limitations (M10).
- End-to-end test fixture in `tests/test_data/e2e/` (11 real 16S sequences, 600 simulated ONT-like reads in known proportions plus 30 non-16S reads) and an `e2e` nf-test that runs the real tools with Docker; CI runs it on pushes to `main`, `dev` and `roadmap/**` and on pull requests into `main` or `dev`.
- Parameters for abundance estimation, which were fixed in the RUN_VI script: `--estimator` (`vb`, the Bayesian VI, by default; `em` for Emu's EM loop, kept as a baseline), `--min_abundance` (0.001, the value used so far), `--vi_tolerance` (1e-7), `--vi_max_iterations` (1,000 CAVI iterations per pruning round), `--alpha_prior` (1.0) and `--vi_elbo_threshold` (1e-4, EM only). The README explains the detection floor: species below max(1/R, min_abundance), or max(10/R, min_abundance) above 1,000 reads, are pruned. With the defaults the results are unchanged (F20, M7).
- `--credible_intervals` (off by default, VB only) adds `ci_lower` and `ci_upper` (95% interval of each species' mean-field Beta marginal) to the per-sample tables; `shared_read_fraction`, which it also added at first, is now always written (see Changed) (F20, D7, D27).
- `--max_cpus` (16), `--max_memory` (128.GB) and `--max_time` (240.h) cap what any task asks for, on every retry (`process.resourceLimits`), so a retry no longer asks for more than a node has. The trace file also records each task's requested `cpus`, `memory` and `time` (F28).
- The minimap2 index is built once and kept in `--index_cache_dir` (default `<output_dir>/index_cache`) under a name that records the database (the first 12 digits of the FASTA's MD5), `--type`, `--kmer_size` and the minimap2 version; a later run with the same settings reuses it instead of rebuilding it. `--db_index` takes a prebuilt `.mmi` and skips the indexing step. `bin/build_gtdb_db.py` no longer suggests indexing with `minimap2 -k 27`. INDEX_DB builds with the CPUs it asks for (`-t`; the index is the same) and stops if the minimap2 it runs is not the version in its cache name. The README no longer tells users to pass a `.mmi` file as `--db`, which the pipeline rejects (F26).
- `--min_quality` (off by default) turns on fastplong's quality filter on mean read quality; without it reads are filtered by length only, as before (F20, D5).
- Parameters are validated against `nextflow_schema.json` (now JSON Schema draft 2020-12, with `$defs`) by the nf-schema plugin (2.8.0) at the start of every run: types, allowed values, ranges and that input paths exist. `main.nf` declares each parameter's type, so values given on the command line, which Nextflow 26.04 passes as strings, are converted: `--keep_counts false` and `--keep_files false` now turn those options off (before, the string "false" counted as true), `--min_length 400` no longer stops the run with a type error, and `--N abc` is rejected (F40).
- Methods benchmark: the scripts in `benchmarks/` and the results in `docs/benchmarks/` (`methods_2026.md`, `m10_claims.md`, tables and figures with their source data) on simulated and real mock communities, against Emu, EPI2ME wf-16s, NanoCLUST and Kraken2 + Bracken. Experimental options used by the benchmark, not exposed as pipeline parameters and off by default: `read_likelihood.py --error-profile`, `--score per-base`, `--min-identity`, `--min-aligned-length` and `--unknown-delta` (per-read floors of the unknown-taxon component), and `estimate_abundance.py --no-pruning` and `--unknown-component`. Pipeline outputs are unchanged (M1-M9, F18). The options that lost were removed after the benchmark; its experimental arms run with `bin/` at commit 34d0411: the unknown-taxon component (`--unknown-delta`, `--unknown-component` and the `unknown_floor` array of the scores file; D22, M5); the identity and aligned-length filters (`--min-identity`, `--min-aligned-length`; D25); the per-base read score (`--score per-base`; Emu's score stays, D29); and `--error-profile all` and `--no-pruning` (the error profile stays primary-only and EM stays as the baseline estimator, D31).
- `bin/build_gtdb_db.py` and `--cmd build-database` write two reports into the database directory, from the new `bin/identical_16s.py` (which also runs on an existing database): `identical_16s.tsv`, the species that share an identical 16S sequence with another species of their genus, with their groups and the species each shares a sequence with (GTDB r226: 4,335 species in 1,264 groups when built with `--merge-gtdb-suffixes`, 5,215 in 1,513 with the default, which keeps suffixed siblings apart; the Phase 7a audit's 3,925 left out the sequences that a species of another genus also carries; D24), and `cross_genus_identical.tsv`, the sequences found in more than one genus, likely contaminant contigs (r226: 1,360 sequences carried by 2,080 species; D23). Genera are compared with GTDB's clade suffixes stripped (one or two capitals after an underscore that end a word, `_A` to `_ZZ`, so other databases' names such as `Lachnospiraceae_NK4A136_group` or `Family_XIII` stay whole). The database itself is unchanged (D23, D24, M6).
- Per-sample tables have an `identical_16s_species` column: for each reported species, the species of the same genus that share an identical 16S sequence with it in the database, which 16S cannot separate (on GTDB r226, 4,335 species, or 5,215 with suffixes kept; for example *Staphylococcus aureus* lists *S. argenteus*, which the Phase 7a simulations estimated at 0 against a true 50% of its pair). It comes from the database's `identical_16s.tsv`; for a database without one, a new step, IDENTICAL_16S, builds it once from the FASTA and taxonomy and keeps it in `--index_cache_dir` (D24, M6, W13).

### Changed

- The README and `docs/output.md` describe what the code does: every parameter with its default (`--keep_counts` is on by default; `--db` is a database directory, never a `.mmi` file), the output folders (`sample_outputs/`, `consolidated_output/`, `qc/fastplong/`, `filtered/`, `bam/`, `index_cache/`), every column of the per-sample table, the four consolidated tables and the database files; the workflow diagram names the length filter, which is what FILTER_READS applies unless `--min_quality` is given (F42, F22).
- NanoVI needs Nextflow 26.04.0 or later (was 24.0.0): the typed parameter declarations and nf-schema 2.8.0 need it. `--K` is a string, so it also takes minimap2's suffixes (`500M`) (F40).
- The samplesheet is read with nf-schema's `samplesheetToList` against `assets/schema_input.json`, replacing `lib/SamplesheetParser.groovy`. The same rules apply (unique sample IDs of letters, digits, `_` and `-`; an existing `.fastq`, `.fq`, `.fastq.gz` or `.fq.gz` file), but FASTQ paths may now be relative to the launch directory or remote (`s3://`, `gs://`, `az://`, `https://`), and quoted fields are read correctly. Errors list every bad row. Lines starting with `#` are no longer skipped: each is read as a sample and fails validation (F15).
- Image tags of the form `X.Y.Z` are published only from `vX.Y.Z` git tags and never overwritten: the image job fails if the tag does not match the top release in `CHANGELOG.md` or if the version already exists on GHCR, and `:latest` moves only with a release. Pushes to `main` publish `:main` and `:sha-<commit>` instead of re-pushing the CHANGELOG version and `:latest` on every merge. A `CHANGELOG.md` without a release heading fails the job instead of falling back to 1.0.0. The image is published only after the version check and every test job (unit, stub, e2e with Docker and with conda) pass. The README has a Releasing section (F34).
- The Python image (`containers/Dockerfile`) is built on Python 3.11 (Python 3.8 is past end of life) with exact versions: numpy 2.4.6, pandas 3.0.6, scipy 1.17.1, pysam 0.24.1, biopython 1.88. `containers/Singularity.def` installs from the same `requirements.txt` and no longer installs flatten-dict. The unit tests run on Python 3.11 and 3.13 (F32).
- Every Python step runs in that image, `ghcr.io/microbialds/nanovi-python:1.0.0`, the one CI builds, instead of the personal `ccuriqueo/vi-python:3.8` (numpy 1.21.6, pandas 1.3.5, no scipy). The e2e CI job builds it from the commit under test. Abundances are unchanged (F30).
- fastplong runs from its BioContainers image instead of `ccuriqueo/fastplong:0.2.2`, and minimap2's image is pinned to a build (versions below). With `-profile singularity` both are downloaded as the matching builds from depot.galaxyproject.org instead of being converted from Docker images (F30).
- The `docker` and `standard` profiles no longer mount the host's `/usr/bin/ps` into the containers. Every image has its own `ps`, and the host binary could not run inside them (it needs the host's `libprocps`), so Nextflow recorded almost no memory use: the trace now shows, for example, 81 MB peak RSS for RUN_VI on the e2e fixture instead of 3 MB (F30).
- `-profile test` now runs on the bundled e2e fixture instead of a `db/` directory that was not in the repository, and uses the pipeline's default k-mer size (21) instead of 27.
- CI runs on pushes to `dev` and `roadmap/**` and on pull requests into `dev`.
- README credits Emu (Curry et al., Nature Methods 2022) for the CIGAR-based read score, the EM loop and the database layout; Emu's MIT license is in `THIRD_PARTY_LICENSES`, and the adapted files name their origin.
- The error profile (the share of insertions, deletions, soft clips and mismatches that weights each read's score) is estimated from primary mapped alignments only, as Emu does; it also counted secondary and supplementary alignments, which go to less similar references and raised the mismatch share. Abundances move slightly (F03).
- The Bayesian VI's inner loop stops when the largest change in any species' posterior mean abundance is below 1e-7 (at most 1,000 iterations) instead of on a relative ELBO change of 5e-4, which stopped large samples after a couple of iterations, before minor species had settled. The ELBO is still computed and a warning is printed if it falls. Abundances move slightly (F05).
- The Bayesian VI's pruning loop stops only when no species is below the pruning threshold. It no longer stops on an ELBO comparison between different species sets, which could leave species below the threshold in the output, and each round starts from the previous solution instead of from uniform. After the last allowed round (10) it fits the pruned set once more and warns if species are still below the threshold (F07).
- The Python steps GET_CIGAR, COMPUTE_LOGP, RUN_VI and WRITE_OUTPUT run scripts in `bin/` (now `read_likelihood.py`, `estimate_abundance.py` and `write_abundance_table.py`, each with `--help`; see F24) instead of inline `python3 -c` programs, and every process calls its script by name from the PATH Nextflow sets up, instead of staging `bin/` as an input. Outputs are unchanged (F38).
- The minimap2 index is built with the alignment preset (`-x <type> -k <kmer_size>`) instead of `-k` alone. Index-time settings (k, window, homopolymer compression) are fixed in the index, so `--type map-pb` aligned with an index built without its homopolymer compression. With the default `--type map-ont` the index is byte-identical to before; results change only for other `--type` values (F19).
- The Python steps ask for one CPU each instead of 2 to 4, so more samples run at once on the same allocation. They also run on one core now: NumPy's OpenBLAS started a thread per core of the node, so a step used several cores (up to 8 on an 80-core machine) whatever it asked for; `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` are set to 1 (`env` in `nextflow.config`), with identical results. fastplong now runs with the CPUs FILTER_READS asks for (`--thread`, 2; it used its default 3 worker threads), with identical output. Alignment keeps 8 CPUs on every retry, and only its memory and time grow: the retried exit codes are memory and time kills, which more threads don't fix (8, 16 and 24 CPUs before) (F28).
- `--K` defaults to `500M` (minimap2's default: 500 million read bases per batch) instead of `4000000000`, which loaded up to 4 billion bases, usually every read of a sample, at once. Alignments do not depend on the batch size (F18).
- fastplong writes the filtered reads gzipped instead of a plain FASTQ that was then compressed in a second pass, and minimap2's output goes through `samtools view -b` into a BAM file instead of a plain SAM (325 KB instead of 2.0 MB on the e2e fixture; LIKELIHOOD reads either). fastplong compresses at its own default level, so the filtered FASTQ is a little larger than before (256 KB instead of 231 KB on the fixture) but is no longer compressed a second time. With `--keep_files` the run publishes `bam/<sample>_vi_alignments.bam` instead of `sam/<sample>_vi_alignments.sam`. ALIGN_READS runs a BioContainers image with both minimap2 and samtools (versions below) and reports both. Every process runs bash with `pipefail`, so a failing step in a pipe fails the task; when minimap2 fails, ALIGN_READS exits with minimap2's status, so a memory kill is still retried with more memory, and when samtools fails first it keeps samtools' status. samtools compresses with two extra threads (F27).
- The Bayesian VI keeps each read's few candidate species in padded reads × candidates arrays instead of dense reads × species matrices, so RUN_VI's memory grows with the number of reads and no longer with reads × species (3.2 GB per array at 200,000 reads and 2,000 species before). Abundances agree with the dense version within 1e-9; on the e2e fixture the tables are identical (F23).
- GET_CIGAR and COMPUTE_LOGP are merged into one step, LIKELIHOOD (`bin/read_likelihood.py`), which reads the alignments once instead of twice. It hands the per-read scores to RUN_VI as a compressed NumPy file (`<sample>_read_scores.npz`) instead of a JSON of per-read dictionaries, and the read counts to WRITE_OUTPUT as a small JSON (`<sample>_read_counts.json`), so WRITE_OUTPUT no longer stages the per-read data. Outputs are unchanged (F24).
- fastplong 0.7.0 (was 0.2.2) and minimap2 2.31 (was 2.24, now with samtools 1.24 instead of 1.18), from BioContainers (`quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0` and the mulled minimap2 + samtools image `mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:ad4c240457048f1d3cc0bd8c8fa39bf27d222eb7-0`), and from bioconda under `-profile conda`. On the Phase 7a real mock communities the upgrade gives the same species, species L1 within 0.0004 and assigned reads within 3 per 10,000 reads. The cached index name carries the minimap2 version, so existing indexes are rebuilt once; an index built with 2.24 and passed as `--db_index` should be rebuilt with 2.31 (D30, D16).
- `bin/build_gtdb_db.py` keeps GTDB's names as GTDB writes them, clade suffixes included (`Listeria monocytogenes_B`, `Pseudomonas_E aeruginosa`), so species and genera that GTDB keeps apart get their own tax IDs. Before, the suffixes were always stripped: on GTDB r226 that merged 5,037 GTDB species into 1,614 names (4,787 into 1,529 among species with a sequence of at least 900 bp), 35 of them across GTDB genera, each merged species keeping the lineage of its first genome. `--merge-gtdb-suffixes` gives the old database (byte-identical on r226) and writes `merged_species.tsv`, which lists every name that merges two or more GTDB species with a sequence in the database (r226: 1,529 names, 26 across genera), its GTDB species and genera, and the ranks at which their lineages differ. A database built with the new default has different tax IDs and more species than one built before, and its tables can report GTDB sibling species separately: in the Phase 7b confirmation *L. monocytogenes* and *L. monocytogenes_B* both appear in 70 of 93 samples, and *Streptococcus mitis* splits in 36 (`docs/benchmarks/methods_2026.md` §13) (F16, D4).
- Species names in the output tables are kept whole, as the database's `taxonomy.tsv` gives them, instead of being cut to their first two words, which merged distinct species with longer names (`Candidatus Pelagibacter ubique` became `Candidatus Pelagibacter`) and NCBI subspecies and strains. `--binomial_names` cuts names to genus and epithet as before, but keeps `Candidatus` names to three words; a row that then sums several tax IDs counts the reads split among them in its `shared_read_fraction`. GTDB names have two words, so GTDB-based tables are unchanged; tables from NCBI or custom databases with longer names change (F17, D4).
- Every per-sample table has a `shared_read_fraction` column, with either estimator: the share of each species' expected reads that it splits with another reported species, a per-species flag of ambiguity (in the Phase 7a benchmark it flagged the estimates that missed the truth, AUC 0.75-0.84). RUN_VI always writes the per-species read summary (`read_summary` in `<sample>_abundance.json`, which replaces `posterior`; EM's comes from its final responsibilities over the species that survived pruning, as VB's comes from its final fit, so the column means the same with either estimator). `--credible_intervals` stays off by default and now adds only `ci_lower` and `ci_upper`, in WRITE_OUTPUT; they are described as sampling intervals, not calibrated uncertainty: their 95% coverage was 93%, 81% and 38% at 1,000, 10,000 and 50,000 reads, and mean-field, a read bootstrap and the exact posterior agreed (D27, D7, M2). Abundances are unchanged.

### Removed

- The standalone `abundance` subcommand of `bin/main.py`, with `bin/alignment.py` and `bin/utils.py` (decision D1). It could not produce results (its tax IDs were strings, so no read matched a species), used its own defaults (k-mer size 27, fastplong's quality filter on) and wrote the minimap2 index into the database directory. Abundance estimation runs only through the Nextflow pipeline; `bin/main.py` keeps `build-database`, `collapse-taxonomy` and `combine-outputs` (F02, F39).
- `output_read_assignments`, which nothing called and which would have failed (it used pandas without importing it) (F09).

### Fixed

- `--keep_files` never published anything: the `enabled` option of the SAM and filtered-FASTQ publish rules was a closure, which Nextflow never evaluates as true. With `--keep_files` the run now writes `bam/` (the alignments, see F27) and `filtered/` (F41).
- `build-database` published the database into `<output_dir>/<db_name>/<db_name>/`; it now goes to `<output_dir>/<db_name>/`. With `--ncbi_taxonomy`, `taxonomy.tsv` lists its rows sorted by tax ID instead of in a different order on each run, so the same inputs give the same file (F41).
- `-profile conda` works: every process declares its conda packages (fastplong, minimap2 and samtools from bioconda, the same builds as the containers), and the Python steps use `containers/environment.yml`, which pins the image's Python and package versions (it replaces the top-level `environment.yml`, which pinned Python 3.8 and lacked fastplong and scipy). Before, no process declared a conda environment, so the profile ran every tool from the host PATH. On the e2e fixture the conda run gives the same tables as Docker, and CI runs the e2e test with `-profile conda` (through micromamba) as well as with Docker (F31).
- `RUN_VI` no longer runs `pip install scipy` in every task: scipy is in the image, so the step works on compute nodes without internet access and every run uses the same scipy. On the e2e fixture RUN_VI takes 2.5 s instead of 17.7 s. `build-database`, `collapse-taxonomy` and `combine-outputs`, which failed under Docker because their image had no scipy, now run (F29).
- The log banner showed `v1.0` and `bin/main.py --version` said 1.0; both now read the pipeline version (1.0.0), from `manifest.version` and from `bin/_version.py`, the Python code's only version string. A CI job fails when `manifest.version`, the top release in this file, the Dockerfile and Singularity.def labels, `bin/_version.py` and the image tag in `conf/base.config` disagree (F33).
- Per-sample tables now keep an `Unassigned` row (abundance 0, `estimated_counts` = reads with no alignment), so `otu_table_counts.tsv` and `sample_metadata.tsv` report unassigned reads instead of 0. Tax IDs missing from `taxonomy.tsv` are reported as `Unknown taxid <id>` with a warning, and species with an empty rank are no longer dropped. The blank ranks of these rows carry their label, so `collapse-taxonomy` and `combine-outputs` keep them apart at every rank (F01).
- Reads whose candidate species were all pruned by the VI are no longer spread over the reported species. Per-sample tables have a `Below detection` row (abundance 0, `estimated_counts` = those reads) just before `Unassigned`, and the reported species share the remaining assigned reads. The consolidated tables carry the row too (blank taxonomy, like `Unassigned`); `assigned_counts` in `sample_metadata.tsv` still counts every read with an alignment. The VI's uniform start counts only reads with a score (F06).
- `estimated_counts` uses largest-remainder rounding instead of truncation, so the assigned rows sum exactly to the number of assigned reads (F08).
- The consolidated tables sum rows that share a species name instead of keeping the last one, take each species' taxonomy from its most common lineage across samples, and skip `otu_table_counts.tsv` with a warning when no sample has counts (`--keep_counts false`) instead of writing a table of zeros (F11).
- `combine-outputs` stops with an error naming the file when a table lacks the rank or metric column (for example `--counts` on tables written without counts), writes 0 instead of a blank for taxa absent from a sample (whole numbers with `--counts`), orders sample columns by file name, and reads only `*_rel-abundance.tsv`, not the tables `collapse-taxonomy` writes next to them (F12).
- `build-database` stops early unless it gets `--sequences`, `--seq2tax` and exactly one of `--taxonomy_list` or `--ncbi_taxonomy`, and `--ncbi_taxonomy` now reaches the database builder (it was declared but never passed). `collapse-taxonomy` and `combine-outputs` check that `--rank` is given and is a valid rank (F13).
- Directory input removes only the FASTQ extension to make sample IDs (`S1.run1.fastq.gz` was `S1`, now `S1.run1`), and stops with a list of the files when two of them give the same ID instead of merging their results, or when an ID has characters other than letters, digits, `.`, `_` and `-` (F14).
- The CIGAR statistics step reads the SAM once in a single streaming pass instead of handing every alignment to a multiprocessing pool that only added up lists (F25). It no longer crashes when a sample has fewer alignments than CPUs (including an empty SAM), keeps each read's longest alignment over all of its alignments instead of the value from the last chunk that saw the read, and stops with a clear error when the SAM has no mapped alignment (F04). Its CPU count no longer matters.
- An alignment that uses an operation absent from the error profile (now possible, since the profile comes from primary alignments only) is skipped as before, but its read is counted as unassigned when it has no other scored alignment instead of disappearing from the totals (F10).
- A sample with no read aligned to the database, such as a blank, or a sample whose reads all fail the length filter, no longer stops the run (LIKELIHOOD failed and no sample got its tables). Its table has only the `Below detection` (0) and `Unassigned` rows, with every read unassigned. The consolidated tables keep it with zero counts, and the run log warns, naming the sample. When no read passes the length filter, `filtered/<sample>_filtered.fastq.gz` is now an empty gzip file rather than fastplong's 0-byte file (F45).
- `shared_read_fraction` is computed per row of the table. A row that sums several tax IDs (with `--binomial_names`, or two tax IDs with the same species name) no longer counts the reads split only among those tax IDs as shared. With `--binomial_names`, two strains that split their reads 50/50 gave the merged species a fraction of 1.0; it is now 0 (0.2 in the test, from the reads it really shares with a sibling). RUN_VI now reads `taxonomy.tsv` and `--binomial_names` to find the rows. Tables without summed rows are unchanged (F46).
- With `--estimator em`, the abundances after pruning sum to 1. The EM inner loop returned the frequencies from before its last update, which do not sum to 1 right after pruning (by up to its 5e-4 tolerance). On the e2e fixture *E. coli* and *S. flexneri* move by 2.6e-4 and the estimates are closer to Emu 3.6.2's on the same scores (L1 2.9e-5 instead of 5.5e-4) (F47). `estimated_counts` divides the abundances by their sum before rounding and hands out every leftover read, so the species rows always add up to the assigned reads; abundances summing to 0.9995 lost about 50 of 100,000 reads, and a sum above 1 over-counted. Counts on the e2e fixture are unchanged (F48).

## [1.0.0] - 2026-06-27

### Added

- **Phyloseq-compatible consolidated output**: The `abundance` command now automatically generates consolidated tables suitable for direct import into R phyloseq or Python analysis
  - New `bin/generate_phyloseq_tables.py` script that reads per-sample `*_rel-abundance.tsv` files and produces:
    - `otu_table_abundance.tsv` — species × samples relative abundance matrix
    - `otu_table_counts.tsv` — species × samples estimated counts matrix
    - `tax_table.tsv` — species × taxonomic ranks (genus, family, order, class, phylum, superkingdom)
    - `sample_metadata.tsv` — per-sample summary (total abundance, taxa detected, assigned/unassigned counts)
  - New Nextflow module `GENERATE_PHYLOSEQ_TABLES` in `modules/local/generate_phyloseq_tables.nf`
  - Wired into the `ABUNDANCE` subworkflow as Step 8, collecting all per-sample outputs before consolidation
  - Includes "Unassigned" row in OTU tables; uses unified 2-word species names
  - Inspired by BugBuster's `TAXONOMY_PHYLOSEQ` module pattern

- **Container build definitions**: Added build files for the NanoVI Python environment so users can rebuild containers and HPC users can create Singularity images
  - Populated `containers/Dockerfile` (FROM python:3.8-slim, system deps via apt, Python deps via pip)
  - Created `containers/requirements.txt` with versioned ranges for numpy, pandas, pysam, biopython, and flatten-dict
  - Populated `containers/Singularity.def` (Bootstrap from same Docker base for consistency, %post, %environment, %runscript)
  - Created `containers/.dockerignore` to exclude pipeline files from the Docker build context
  - Added "Container Images" section to `README.md` with Docker and Singularity build commands and custom config override instructions
  - Note: FastpLong (`ccuriqueo/fastplong:0.2.2`) and minimap2 (`quay.io/biocontainers/minimap2:2.24--h7132678_1`) are pulled from public registries and do not require custom build definitions

- **Parameter schema and runtime validation**: Added `nextflow_schema.json` and hardened `main.nf` abundance block with early-exit checks
  - Created complete `nextflow_schema.json` with four definition groups (`general_options`, `abundance_options`, `database_options`, `postprocessing_options`), covering all pipeline parameters with types, defaults, enums, and min/max constraints
  - Enables `nextflow run main.nf --help` to display formatted parameter documentation
  - Added null guards for `--input` and `--db` in the abundance block (previously caused cryptic NullPointerException tracebacks)
  - Added `taxonomy_tsv` fallback: defaults to `${params.db}/taxonomy.tsv` when not explicitly provided
  - Added `min_length < max_length` cross-parameter validation with a descriptive error message

- **Samplesheet-based input support**: Implemented nf-core-style CSV samplesheet input for improved metadata handling and validation
  - Created `lib/SamplesheetParser.groovy` for robust CSV parsing with comprehensive validation
  - Added `assets/schema_input.json` JSON Schema for samplesheet validation
  - Updated `assets/samplesheet.csv` with example format
  - Modified `main.nf` to support both samplesheet (`.csv`) and directory input modes
  - Updated `FILTER_READS` process to accept `tuple(sample_id, fastq)` input instead of deriving sample IDs from filenames
  - Validation features:
    - Enforces unique sample IDs (duplicates are rejected)
    - Requires absolute paths for FASTQ files
    - Validates sample names (alphanumeric, underscores, hyphens only)
    - Checks file existence and valid FASTQ extensions
  - Backward compatible: Directory-based input continues to work for quick tests
  - Benefits: Explicit metadata, early error detection, flexible file naming, nf-core compliance
  - Migration: Users can now provide `--input samplesheet.csv` instead of `--input /path/to/fastq/`

### Changed

- **Restructured output directory layout**: Per-sample and consolidated results are now organized into separate folders
  - `WRITE_OUTPUT` publishDir moved from `${params.output_dir}/` to `${params.output_dir}/sample_outputs/`
  - New `GENERATE_PHYLOSEQ_TABLES` publishDir at `${params.output_dir}/consolidated_output/`
  - QC, filtered, SAM, and pipeline_info directories remain unchanged

- **`estimated_counts` column rename, fix, and default enablement**: Corrected the output abundance table column for estimated read counts
  - Renamed column `"estimated counts"` → `"estimated_counts"` in `bin/taxonomy.py` (`freq_to_lineage_df`, `collapse_rank`, `combine_outputs`) for consistency with snake_case naming conventions
  - Fixed empty-column bug in `freq_to_lineage_df`: the column was unconditionally added to `column_order` as an empty string `""` even when `counts=False`; it is now only included when `counts=True`
  - Changed `keep_counts` default from `false` to `true` in `nextflow.config` so estimated counts are included in every pipeline run by default

- **Code quality cleanup across Python scripts and Nextflow modules**
  - `bin/variational_inference.py`: Removed 12 unreachable `pass` statements after `return`, removed unused `from Bio.Seq import Seq` import, added docstrings to 6 public functions (`get_cigar_op_log_probabilities`, `log_prob_rgs_dict`, `variational_inference`, `variational_inference_iterations`, `update_q_distribution`, `maximize_elbo`), improved `ValueError` messages in `variational_inference_iterations` to include actionable guidance
  - `bin/main.py`: Commented out dead `output_read_assignments` call (the `--keep-read-assignments` flag has no effect since `variational_inference_iterations` always returns `None` for `read_dist`)
  - `bin/alignment.py`, `bin/utils.py`, `bin/taxonomy.py`: Removed emoji characters (`🔄`, `✅`, `🧬`) from all output messages
  - `modules/local/align_reads.nf`: Translated Spanish echo message to English; removed emoji
  - `modules/local/filter_reads.nf`: Translated Spanish comment to English

- **Idiomatic `bin/` path resolution**: Replaced fragile `cp -r ${projectDir}/bin bin` shell copies and `sys.path.insert(0, "./bin")` assumptions with Nextflow-native patterns
  - `BUILD_DB`, `COLLAPSE`, `COMBINE`: added formal `path bin_dir` input; `bin_ch` is created once in `main.nf` via `channel.value(file("${projectDir}/bin"))` and threaded through `BUILD_DATABASE`, `COLLAPSE_TAXONOMY`, and `COMBINE_OUTPUTS` subworkflows
  - `GET_CIGAR`, `COMPUTE_LOGP`, `RUN_VI`, `WRITE_OUTPUT`: replaced `sys.path.insert(0, "./bin")` with `sys.path.insert(0, "${projectDir}/bin")` (Nextflow interpolates `projectDir` in script blocks)
  - Removed stray `nextflow.enable.dsl = 2` declarations from `build_database.nf` and `collapse_taxonomy.nf` module files (only `main.nf` requires this)
  - Added missing `output:` block to `BUILD_DB` process (`path "${params.db_name}/**"`)
  - Updated all remaining `Channel.value(...)` calls in `main.nf` to lowercase `channel.value(...)` per Nextflow DSL2 conventions
  - Impact: Fixes Nextflow caching (no undeclared `projectDir` side-effects in process scripts), enables `nextflow run microbialds/NanoVI` remote pull execution, and removes reliance on working-directory layout assumptions

- **Configuration overhaul with multi-profile support and error handling**: Complete restructuring of pipeline configuration for improved portability, robustness, and maintainability
  - Rewrote `nextflow.config` with comprehensive profile system:
    - Added `docker`, `singularity`, `conda` profiles for container engine selection
    - Added `slurm`, `sge` profiles for HPC execution
    - Added `test` profile with minimal resources for CI/testing
    - Added `standard` profile as default (Docker-based) for backward compatibility
  - Rewrote `conf/base.config` with robust error handling and dynamic resource allocation:
    - Implemented automatic retry strategy for common failure exit codes (104, 134, 137, 139, 143, 247)
    - Set `maxRetries = 2` and `maxErrors = '-1'` for resilient execution
    - Converted all resource allocations to use dynamic scaling with `task.attempt` multiplier
    - Updated memory/time syntax to Nextflow standard format (e.g., `2.GB`, `30.min`)
    - Added dynamic CPU scaling for alignment process
  - Created `conf/modules.config` for centralized output publishing:
    - Moved all `publishDir` directives from individual process files to centralized configuration
    - Configured publishing for 6 processes: FILTER_READS, ALIGN_READS, WRITE_OUTPUT, BUILD_DB, COLLAPSE, COMBINE
    - Implemented conditional publishing based on `params.keep_files` flag
  - Created `conf/test.config` with reduced resource allocations for testing
  - Created `environment.yml` for conda profile support with all required dependencies
  - Added pipeline manifest with metadata (name, author, description, version, homepage)
  - Enabled automatic generation of timeline, report, and trace files in `${params.output_dir}/pipeline_info/`
  - Impact: Pipeline now supports multiple execution environments (local, HPC, cloud), automatically retries failed tasks, scales resources on retry, and provides centralized configuration management
  - **Note**: This supersedes previous configuration-related changes from prompts 05 (portability), 07 (error handling), and parts of 10 (validation)

- **Refactored to DSL2 subworkflows architecture**: Extracted monolithic workflow logic into composable, testable named subworkflows
  - Created `subworkflows/local/abundance.nf` with complete 7-step abundance pipeline (filter → align → CIGAR → log-prob → VI → output)
  - Created `subworkflows/local/build_database.nf` as wrapper for database building process
  - Created `subworkflows/local/postprocessing.nf` with two workflows: `COLLAPSE_TAXONOMY` and `COMBINE_OUTPUTS`
  - Rewrote `main.nf` as clean dispatcher (~95 lines vs ~109 lines) that includes and calls subworkflows
  - Moved module includes from top-level to subworkflow-level for better encapsulation
  - Added input validation for abundance subcommand (database directory, FASTA file, taxonomy file existence checks)
  - Added version logging banner and improved error messages
  - Impact: Improved code organization, testability, and maintainability while maintaining identical functionality and CLI interface
  - Each subworkflow can now be tested independently and the codebase follows Nextflow DSL2 best practices

- **CPU allocation consistency and parameter rename**: Standardized CPU resource management across all processes
  - Renamed pipeline parameter from `--threads` to `--cpus` in `nextflow.config` for consistency with Nextflow's resource allocation model
  - Modified `GET_CIGAR` process in `modules/abundance/cigar_probs.nf` to use `${task.cpus}` instead of accepting `val threads` parameter
  - Updated `main.nf` to remove `params.threads` argument from `GET_CIGAR` call (line 42)
  - Adjusted CPU allocations in `config/containers.config`:
    - `alignment` label: 10 → 8 CPUs (aligned with default params.cpus)
    - `vi_stats` label: 2 → 4 CPUs (compromise for GET_CIGAR multiprocessing and COMPUTE_LOGP single-thread)
  - Updated `README.md` documentation to reflect new `--cpus` parameter name
  - Impact: Ensures tool parallelism matches Nextflow scheduler resource allocation, preventing CPU oversubscription/underutilization
  - **Breaking change**: Users must now use `--cpus` instead of `--threads` when running the pipeline

- **Performance optimization in minimap2 indexing**: Refactored alignment workflow to build database index once and reuse across all samples
  - Created new `INDEX_DB` process in `modules/abundance/index_db.nf` to build minimap2 index as a separate step
  - Modified `ALIGN_READS` process in `modules/abundance/alignment.nf` to accept pre-built index instead of building it internally
  - Removed conditional indexing logic (lines 16-19) that was always rebuilding the index due to fresh Nextflow work directories
  - Updated `main.nf` to wire `INDEX_DB` before `ALIGN_READS`, ensuring index is built once and cached
  - Also fixed thread allocation in `ALIGN_READS`: changed from `${params.cpus}` to `${task.cpus}` for proper resource management
  - Impact: Significant time savings when processing multiple samples, especially with large databases like GTDB (indexing now happens once instead of per-sample)

- **Parameter default reconciliation**: Aligned all parameter defaults to a single source of truth across `nextflow.config`, `bin/main.py`, and `README.md`
  - `bin/main.py`: Fixed `--kmer-size` default from `28` to `27` (matches `nextflow.config`)
  - `bin/main.py`: Fixed `lli_thresh` (VI convergence threshold) from `0.01` to `0.0001` (matches `RUN_VI` module)
  - `README.md`: Fixed `--N` default from `20` to `3` (matches `nextflow.config` and `main.py`)
  - `README.md`: Fixed `--K` default from `1000000000` to `4000000000` (matches `nextflow.config` and `main.py`)
  - `README.md`: Added deprecation note for `--threads` in parameter table (removed as user-facing parameter in prompt 04; replaced by `--cpus`)
  - `nextflow.config`: Verified clean — `params.rank` appears exactly once; all values already correct

### Fixed

- **Removed `flatten_dict` dependency from `bin/variational_inference.py`**: The container `ccuriqueo/vi-python:3.8` did not have `flatten-dict` installed, causing `GET_CIGAR` to fail silently with exit status 1
  - Replaced `from flatten_dict import unflatten` import and its usage in `update_q_distribution` with an equivalent manual nested dict construction
  - Removed `flatten-dict` from `containers/requirements.txt` and `environment.yml`

- **Critical bug in `bin/variational_inference.py`**: Fixed incorrect dictionary key check in `process_alignment_chunk` function (line 38)
  - Changed condition from `if align_len not in dict_longest_align:` to `if query_name not in dict_longest_align:`
  - Previously, the dictionary was incorrectly indexed by alignment lengths (integers) instead of read names (strings)
  - This caused incorrect length-normalization of log-likelihoods in downstream `compute_log_prob_rgs` function
  - Impact: Ensures proper tracking of longest alignment per read for accurate taxonomic classification

- **Race condition in `WRITE_OUTPUT` process**: Fixed data mismatch when processing multiple samples concurrently
  - Modified `main.nf` (lines 47-67) to join `freq_output` and `logp_data` channels by sample_id before calling `WRITE_OUTPUT`
  - Rewrote `modules/abundance/write_output.nf` to accept tuple(sample_id, abundance_json, logp_json) and extract counts internally
  - Previously, `assigned` and `unassigned` counts were passed as separate val channels without sample_id association
  - This caused non-deterministic pairing where counts from sample A could be matched with abundance from sample B
  - Impact: Ensures correct correspondence between abundance data and read counts for each sample
