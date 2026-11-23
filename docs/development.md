# Developing NanoVI

For people changing NanoVI's code or images, or making a release. Running NanoVI needs none
of this: see the [README](../README.md). The repository's layout and conventions (where a new
parameter goes, how modules call the Python scripts) are in `CLAUDE.md` at its root.

## Building the containers

Every Python step runs in one image, `ghcr.io/microbialds/nanovi-python:<version>` (Python 3.11, with the exact package versions in `containers/requirements.txt`); nothing is installed at run time, so the pipeline runs on compute nodes without internet access. To test a change to the image, build it locally under the tag `conf/base.config` uses, and Docker runs the local copy instead of pulling it:

**Docker:**
```bash
docker build -t ghcr.io/microbialds/nanovi-python:1.2.0 -f containers/Dockerfile containers/
```

**Singularity** (from the `containers/` directory, so the build finds `requirements.txt`):
```bash
cd containers && singularity build nanovi-python.sif Singularity.def
```

The run report (`--report`) has its own image, `ghcr.io/microbialds/nanovi-report:<version>`, which installs `containers/report/environment.yml` (conda-forge: Python 3.13, Typst, Jinja, matplotlib, Altair, vl-convert and the Inter fonts) with micromamba; build it the same way:
```bash
docker build -t ghcr.io/microbialds/nanovi-report:1.2.0 containers/report/
cd containers/report && singularity build nanovi-report.sif Singularity.def
```

## Testing

The unit tests cover the Python code the pipeline runs (likelihood, both estimators, output tables, database builds and reports, and the `bin/` scripts); the pipeline tests run `main.nf` with nf-test.

**Unit tests** (pytest, no database needed; the PDF test needs Typst and is skipped without it):
```bash
pip install -r containers/requirements.txt -r tests/requirements-report.txt pytest
pytest tests/unit/ -v
```

**Pipeline tests** (nf-test; each test is tagged `stub` or `e2e`):
```bash
nf-test test tests/main.nf.test --tag stub   # stub mode: channel wiring and parameter errors, no tools or database needed
nf-test test tests/main.nf.test --tag e2e    # real tools with Docker on the fixture in tests/test_data/e2e
NANOVI_E2E_PROFILE=conda nf-test test tests/main.nf.test --tag e2e    # the same test with -profile conda
```

The `e2e` tests run fastplong, minimap2 and the Python steps on a small real 16S database (11 sequences, including the near neighbours *E. coli* / *S. flexneri* and *S. pneumoniae* / *S. mitis*) with 600 simulated reads in known proportions plus 30 non-16S reads, and check the result against `tests/test_data/e2e/expected.tsv` (every present species reported, those at 10% or more within 5 percentage points, absent species below 1%), the `Below detection` and `Unassigned` rows, the four consolidated tables and the tool versions, that samples with no aligned read (only non-16S reads, or every read below `--min_length`) get an all-Unassigned table beside a normal sample, the per-read file, the group and cell-level tables and the run report (its PDF's pages and fonts, no network request in the HTML, identical files on a second render), and that `--cmd build-database --gtdb_ssu` rebuilds the fixture's database exactly from its GTDB-format FASTA; see [tests/test_data/e2e/README.md](../tests/test_data/e2e/README.md). `-profile test,docker` runs the pipeline on the same fixture (`-profile test` on its own sets no container engine).

`tests/check_workflow_diagram.sh` checks that [docs/workflow.md](workflow.md) names every process and subworkflow the code includes, and no other. `docs/images/render_overview.sh` refreshes the embedded fonts and the PNGs of the two figures, `docs/images/nanovi_overview.svg` (the README's) and `docs/images/nanovi_workflow.svg` (docs/workflow.md's), after either SVG is edited; both are hand-written, so a new process or option needs a box in `nanovi_workflow.svg` (and in `nanovi_overview.svg` if a user would notice it).

Continuous integration (GitHub Actions) runs on every push to `main`, `dev`, `feat/**`, `fix/**` and `roadmap/**`, on `v*` tags and on every pull request into `main` or `dev`: the unit tests on Python 3.11 and 3.13, a check that the version strings agree and that docs/workflow.md names every process, the `stub` tests, and the `e2e` tests with Docker (on the Python and report images built from the commit, then the report's tests inside its image, where Typst is) and with conda. Both Docker images are built on every run; they are published to `ghcr.io/microbialds/nanovi-python` and `ghcr.io/microbialds/nanovi-report` only from `main` (as `:main` and `:sha-<commit>`) and from release tags (see Releasing).

## Releasing

Image tags of the form `X.Y.Z` come only from `vX.Y.Z` git tags and are never overwritten: to change a released image, release a new version.

1. On `dev`, rename `## [Unreleased]` in `CHANGELOG.md` to `## [X.Y.Z] - YYYY-MM-DD` (and start a new empty `[Unreleased]` section), and set the same version in `nextflow.config` (`manifest.version`), `containers/Dockerfile` (`LABEL version`), `containers/Singularity.def` (`Version`), `bin/_version.py`, the `nanovi-python` image tags in `conf/base.config`, and for the report image `containers/report/Dockerfile` (`LABEL version`), `containers/report/Singularity.def` (`Version`) and its `nanovi-report` tag in `conf/base.config`. The `Version strings agree` CI job fails if any of them differ. Also update the `-r vX.Y.Z` of the README's examples, the build commands above, and the release and image tag of the report command in `docs/output.md` (Rendering a report from a finished run).
2. Merge `dev` into `main` through a pull request. CI publishes `:main` and `:sha-<commit>`.
3. Tag the merge commit on `main` and push the tag:
   ```bash
   git tag -a vX.Y.Z -m "NanoVI X.Y.Z"
   git push origin vX.Y.Z
   ```
   CI runs the tests and publishes `:X.Y.Z` and `:latest` of both images. It stops without publishing if the tag does not match the top release in `CHANGELOG.md` or if `:X.Y.Z` already exists.
4. Create the GitHub release from the tag, with that version's CHANGELOG section as its notes.

Until step 3 has published `:X.Y.Z`, `-profile docker` and `-profile singularity` cannot pull the images that `conf/base.config` names, so tag soon after merging. A new GHCR package is private: after its first publish, make `nanovi-report` public in its package settings, as `nanovi-python` is.
