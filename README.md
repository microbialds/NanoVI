# **NanoVI**

Pipeline for estimating relative abundances of Oxford Nanopore 16S reads using a variational inference approach. The workflow is written in Nextflow and relies on
Python helper scripts found in the `bin/` directory.


This pipeline provides four main subcommands:

**Abundance:** Full workflow for read filtering, alignment, probability calculation, and abundance estimation.

**Build-database:** Build a custom GTDB-based database from input sequences and taxonomic mappings.

**Collapse-taxonomy:** Collapse a single-sample abundance table to a specified taxonomic rank.

**Combine-outputs:** Merge multiple relative-abundance tables into one multi-sample table at a given rank.

**1. Installation**

Clone this repository, and install:

```
git clone https://github.com/microbialds/NanoVI
cd NanoVI

- [**Nextflow**](https://www.nextflow.io/docs/latest/install.html) v24.x or later
- [**Docker**](https://docs.docker.com/engine/install/)
- Optional environment variable `GTDB_DATABASE_DIR` pointing to a GTDB database
  (used as default for the `--db` parameter)

```

**2. Input file**

NanoVI accepts input files with the extensions .fastq or .fastq.gz. If you need to process multiple samples simultaneously, you should provide the path to the input folder.


**3. Test dependencies**

It looks like you want to test the dependencies of your project using pytest and expect the output to confirm that 9 dependencies have passed.
Navigate to: 

```
cd VI-pipeline 
pytest tests/

```
**4. Command example for Abundance Estimation**




