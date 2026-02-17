# Changelog

All notable changes to the NanoVI pipeline will be documented in this file.

## [Unreleased]

### Changed

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
