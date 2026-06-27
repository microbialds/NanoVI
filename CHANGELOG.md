# Changelog

All notable changes to the NanoVI pipeline will be documented in this file.

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

### Changed

- **Restructured output directory layout**: Per-sample and consolidated results are now organized into separate folders
  - `WRITE_OUTPUT` publishDir moved from `${params.output_dir}/` to `${params.output_dir}/sample_outputs/`
  - New `GENERATE_PHYLOSEQ_TABLES` publishDir at `${params.output_dir}/consolidated_output/`
  - QC, filtered, SAM, and pipeline_info directories remain unchanged

### Changed

- **`estimated_counts` column rename, fix, and default enablement**: Corrected the output abundance table column for estimated read counts
  - Renamed column `"estimated counts"` → `"estimated_counts"` in `bin/taxonomy.py` (`freq_to_lineage_df`, `collapse_rank`, `combine_outputs`) for consistency with snake_case naming conventions
  - Fixed empty-column bug in `freq_to_lineage_df`: the column was unconditionally added to `column_order` as an empty string `""` even when `counts=False`; it is now only included when `counts=True`
  - Changed `keep_counts` default from `false` to `true` in `nextflow.config` so estimated counts are included in every pipeline run by default

### Fixed

- **Removed `flatten_dict` dependency from `bin/variational_inference.py`**: The container `ccuriqueo/vi-python:3.8` did not have `flatten-dict` installed, causing `GET_CIGAR` to fail silently with exit status 1
  - Replaced `from flatten_dict import unflatten` import and its usage in `update_q_distribution` with an equivalent manual nested dict construction
  - Removed `flatten-dict` from `containers/requirements.txt` and `environment.yml`

### Added

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

### Changed

- **Code quality cleanup across Python scripts and Nextflow modules**
  - `bin/variational_inference.py`: Removed 12 unreachable `pass` statements after `return`, removed unused `from Bio.Seq import Seq` import, added docstrings to 6 public functions (`get_cigar_op_log_probabilities`, `log_prob_rgs_dict`, `variational_inference`, `variational_inference_iterations`, `update_q_distribution`, `maximize_elbo`), improved `ValueError` messages in `variational_inference_iterations` to include actionable guidance
  - `bin/main.py`: Commented out dead `output_read_assignments` call (the `--keep-read-assignments` flag has no effect since `variational_inference_iterations` always returns `None` for `read_dist`)
  - `bin/alignment.py`, `bin/utils.py`, `bin/taxonomy.py`: Removed emoji characters (`🔄`, `✅`, `🧬`) from all output messages
  - `modules/local/align_reads.nf`: Translated Spanish echo message to English; removed emoji
  - `modules/local/filter_reads.nf`: Translated Spanish comment to English

### Added

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
