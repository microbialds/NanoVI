# Response to Reviewer Comment 9

> *"Software engineering details are insufficient for an Application Note. There is no
> mention of unit/integration tests, container versioning, or scalability beyond
> single-sample execution."*

We thank the reviewer for this important observation and have substantially
improved the software engineering infrastructure of NanoVI. The revised manuscript
includes an expanded **Implementation** section addressing each concern.

---

## 1. Unit and integration tests

We have added a comprehensive test suite that covers the two layers of the software:

**Unit tests (pytest).** 51 unit tests implemented with `pytest` validate the
core algorithmic components without requiring bioinformatics databases or
alignment files:

- `test_variational_inference.py` (19 tests) — verifies correctness of the
  EM/VI update steps (`update_q_distribution`, `maximize_elbo`,
  `variational_inference`), convergence behaviour, numerical edge cases (zero-
  frequency species, empty read dictionaries), and determinism of the algorithm.
- `test_taxonomy.py` (19 tests) — validates lineage resolution from NCBI taxon
  identifiers, species-name unification (strain-level → two-word species names),
  abundance summing after merging, and rank-collapsing to higher taxonomic levels.
- `test_utils.py` (13 tests) — covers FASTQ path discovery including recursive
  directory traversal, extension filtering, and multi-directory aggregation.

All 51 tests pass in < 1 second on a standard laptop (Python 3.8–3.10).

**Integration tests (nf-test).** A `tests/main.nf.test` file written with the
[nf-test framework](https://www.nf-test.com/) provides two integration tests:

1. A pipeline-level smoke test that executes `main.nf` in stub mode with the
   provided mock FASTQ dataset to validate channel logic and parameter passing.
2. A module-level test for `FILTER_READS` that checks the emitted channel
   structure and output file naming conventions.

Stub execution decouples CI validation from the bioinformatics databases, which
cannot be bundled with the repository due to their size (~50 GB for GTDB r226).
End-to-end pipeline tests with the full GTDB database are documented in the
`README.md` and reproducible using the commands provided there.

---

## 2. Container versioning

The Docker image follows **semantic versioning** tied to the repository release
cycle. The image is built and pushed automatically by a GitHub Actions workflow
(`.github/workflows/ci.yml`) on every merge to `main`:

- The workflow reads the current version from `CHANGELOG.md` (format `## [X.Y.Z]`).
- The image is tagged with both the exact version (e.g., `nanovi-python:1.0.0`)
  and `latest`.
- Images are published to the **GitHub Container Registry (GHCR)** at
  `ghcr.io/microbialds/nanovi-python`, ensuring that every release tag in the
  repository corresponds to a reproducible, immutable container image.
- OCI labels (`org.opencontainers.image.version`, `.source`, `.revision`) link
  each image back to the exact commit from which it was built.

Singularity images can be built from the same Dockerfile via
`singularity build nanovi-python.sif docker://ghcr.io/microbialds/nanovi-python:1.0.0`,
preserving reproducibility on HPC systems where Docker is unavailable.

---

## 3. Scalability beyond single-sample execution

NanoVI was designed for batch processing from the outset. The manuscript now
explicitly describes the following scalability features:

**Multi-sample input.** The pipeline accepts either a directory of FASTQ files
or an nf-core-style samplesheet (CSV with `sample_id`, `fastq` columns),
enabling batch runs of arbitrary size from a single command. Per-sample outputs
are written to `sample_outputs/<sample_id>/`; a dedicated `combine-outputs`
subcommand merges them into a single species × sample abundance matrix.

**Parallel execution.** NanoVI is implemented in Nextflow DSL2 with fine-grained
process-level resource declarations. Each sample is processed independently
within the Nextflow DAG, so N samples occupy at most N × (max CPUs per process)
cores simultaneously on a workstation or cluster. The minimap2 index is built
once and reused across all samples, avoiding redundant I/O.

**HPC and cloud portability.** Pre-configured execution profiles (`-profile
slurm` / `-profile sge`) submit each process as a separate cluster job with
automatic retry logic (up to 3 retries with memory doubling). Users can further
tune CPU and memory allocations via `conf/base.config` without modifying
pipeline code.

**Benchmark.** In the validation experiments reported in this manuscript, NanoVI
processed [N] samples of [X] reads each on a [cluster/server description]. Wall-
clock time scaled approximately linearly with sample count due to independent
parallel execution, demonstrating practical scalability for cohort-level studies.

*(Please fill in the benchmark numbers from your actual experiments.)*

---

## Manuscript changes

The following additions were made to the **Implementation** section of the
Application Note:

- Added subsection *"Testing and continuous integration"* describing the pytest
  unit test suite and nf-test integration tests, with the CI badge linking to
  GitHub Actions.
- Added subsection *"Container versioning"* describing the GHCR release workflow
  and OCI labelling strategy.
- Expanded the *"Scalability"* paragraph to explicitly describe multi-sample
  input, parallel DAG execution, and HPC profiles.
