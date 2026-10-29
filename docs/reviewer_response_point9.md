# Response to Reviewer Comment 9

> *"Software engineering details are insufficient for an Application Note. There is no
> mention of unit/integration tests, container versioning, or scalability beyond
> single-sample execution."*

We thank the reviewer for this important observation and have substantially
improved the software engineering infrastructure of NanoVI. The revised manuscript
includes an expanded **Implementation** section addressing each concern.

---

## 1. Unit and integration tests

The test suite covers the two layers of the software: the Python code that every
pipeline step runs, and the Nextflow pipeline itself.

**Unit tests (pytest).** 227 tests in seven files validate the algorithmic components
without a reference database; they run in under a minute on one core:

| File | Tests | What they check |
|---|---|---|
| `test_variational_inference_bayesian.py` | 53 | The Bayesian estimator (below) |
| `test_em_estimator.py` | 22 | Emu's EM loop, kept as a baseline (`--estimator em`): responsibilities, the frequency update, convergence and pruning, the below-detection count, and agreement with the Bayesian estimator |
| `test_likelihood_io.py` | 19 | The read score on small SAM files written with pysam: the error profile from primary alignments only, the longest alignment over all of a read's alignments, the best score per species, empty or unmapped-only inputs, and the one-pass reader against the two-pass reference |
| `test_outputs.py` | 73 | The per-sample table (Unassigned and Below detection rows, count rounding, species names, the shared-read and interval columns), the consolidated phyloseq tables, `combine-outputs` and the database builders |
| `test_taxonomy.py` | 26 | Lineages from the NCBI taxonomy, species names (whole by default, genus and epithet on request, `Candidatus` names kept whole) and collapsing to higher ranks |
| `test_identical_16s.py` | 8 | The report of species that share an identical 16S sequence and of sequences shared across genera |
| `test_scripts.py` | 26 | Every executable in `bin/`: shebang and permissions, `--help`, and the abundance and table scripts end to end |

**Tests of the Bayesian module.** `test_variational_inference_bayesian.py` tests each
step of the coordinate-ascent variational inference (CAVI) in
`bin/variational_inference_bayesian.py` against an independent calculation:

- *Likelihood storage*: the padded reads × candidates arrays hold the same values as the
  dense reads × species matrix, keep a species' best score when a read aligns to it
  twice, ignore species outside the model and give reads without candidates an empty
  row; pruning turns removed species into padding.
- *Dirichlet expectation*: E[ln π_s] = ψ(α_s) − ψ(Σα) against the digamma formula and a
  Monte Carlo estimate.
- *Responsibilities*: the softmax by hand, rows that sum to one (zero for a read with no
  candidate) and numerical stability for very negative log-likelihoods.
- *Dirichlet update*: α equals the prior plus the expected read counts.
- *Evidence lower bound (ELBO)*: equal to a hand calculation, and never decreasing over
  CAVI iterations.
- *Inner loop*: it stops on the largest change in posterior mean abundance or at the
  iteration cap, stops after one iteration from a converged start, and warns if the
  ELBO falls.
- *Outer loop with pruning*: it recovers known proportions from about 2,000 synthetic
  reads within 0.02, prunes species without support, never stops with species below
  the pruning threshold, warm-starts each round, refits after its last allowed pruning,
  and keeps the most abundant species when every species is below the threshold. The
  result is a fixed point of the updates.
- *Read accounting*: reads whose candidate species were all pruned are counted as below
  detection, not as assigned; the starting point counts only reads with a score.
- *Read summary*: expected and shared reads per species by hand, matching the reported
  abundances, with responsibilities on unreported species left out.
- *Posterior utilities*: Dirichlet entropy against SciPy, Kullback-Leibler divergence
  against a Monte Carlo estimate.
- *Agreement with the dense reference*: every step and the whole estimator agree within
  1e-9 with the earlier dense implementation, kept in `tests/unit/dense_vb_reference.py`,
  including runs with pruning, many species and shuffled species order.

`test_em_estimator.py` also checks that the Bayesian estimator and the EM loop give the
same abundances on synthetic reads (within 0.001 of each other and 0.02 of the truth),
and `test_scripts.py` that both give the same `shared_read_fraction`.

**Integration tests (nf-test).** `tests/main.nf.test` holds 26 tests written with the
[nf-test framework](https://www.nf-test.com/):

1. 24 tests in stub mode, which need no tools or database: the channel wiring and
   published outputs of the whole pipeline, samplesheet and directory input, the
   minimap2 index cache (a second run does not rebuild it), resource caps, and one test
   for each parameter or input error (wrong types, invalid values, duplicate or
   malformed sample IDs, missing files).
2. 2 end-to-end tests that run the real tools (fastplong, minimap2, samtools and the
   Python steps) on a bundled fixture: 11 full-length 16S sequences from distinct
   species, including two near-neighbour pairs, and 600 simulated reads in known
   proportions plus 30 non-16S reads (`tests/test_data/e2e`). They check that the
   assigned abundances sum to 1, that every species in the sample is reported and the
   abundant ones within 5 percentage points of the truth, that absent species stay
   below 1%, that every read is accounted for in the table, the consolidated tables and
   the recorded software versions, and that command-line options such as
   `--keep_counts false` take effect.

**Continuous integration.** GitHub Actions runs the unit tests on Python 3.11 and 3.13,
a check that the version strings agree, the stub tests, and the end-to-end tests twice
(with Docker, on the Python image built from the commit under test, and with conda) on
every push to the main and development branches and on every pull request into them.

---

## 2. Container versioning

Every Python step runs in one image, `ghcr.io/microbialds/nanovi-python:<version>`
(Python 3.11 with exact package versions from `containers/requirements.txt`; nothing is
installed at run time). fastplong, minimap2 and samtools come from pinned BioContainers
images (from bioconda with `-profile conda`).

- Image tags of the form `X.Y.Z` are published to the GitHub Container Registry only
  from `vX.Y.Z` release tags of the repository, and never overwritten: the release job
  stops if the tag does not match the top release in `CHANGELOG.md` or if the version
  already exists. `:latest` moves only with a release; pushes to `main` publish `:main`
  and `:sha-<commit>`.
- A CI job fails when the version in `nextflow.config`, `CHANGELOG.md`, the Dockerfile,
  the Singularity definition, the Python package and the image tag in the pipeline
  configuration disagree.
- OCI labels (`org.opencontainers.image.version`, `.source`, `.revision`) link each
  published image to the commit it was built from.
- Each run records the versions of every tool and Python package it used in
  `pipeline_info/software_versions.yml`.

Singularity images can be built from `containers/Singularity.def`, which installs from
the same requirements file, or, once a release has published the image, pulled from it.

---

## 3. Scalability beyond single-sample execution

NanoVI was designed for batch processing from the outset. The manuscript now
explicitly describes the following scalability features:

**Multi-sample input.** The pipeline accepts either a directory of FASTQ files or an
nf-core-style samplesheet (CSV with `sample` and `fastq` columns; local or remote
files), enabling batch runs of arbitrary size from a single command. Per-sample tables
are written to `sample_outputs/<sample>_rel-abundance.tsv`, and every run also writes
species × sample abundance and count matrices, a taxonomy table and sample metadata in
the layout phyloseq imports (`consolidated_output/`).

**Parallel execution.** NanoVI is implemented in Nextflow DSL2 with process-level
resource declarations. Each sample is processed independently within the Nextflow DAG,
so samples run in parallel up to the resources available. The minimap2 index is built
once and cached under a name that records the database checksum, the alignment preset,
the k-mer size and the minimap2 version, so later runs reuse it. The inference step's
memory grows with the number of reads, not reads × species (200,000 reads over 2,000
species: 58 MB instead of 24 GB).

**HPC and cloud portability.** Pre-configured execution profiles (`-profile
slurm` / `-profile sge`) submit each process as a separate cluster job. A task killed
for lack of memory or time is retried up to twice with more memory and time (the CPU
request stays the same), and `--max_cpus`, `--max_memory` and `--max_time` cap every
request. Users can further tune allocations in `conf/base.config` without modifying
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
