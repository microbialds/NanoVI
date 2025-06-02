# **NanoVI**

Pipeline for estimating relative abundances of Oxford Nanopore 16S reads using a variational inference approach. The workflow is written in Nextflow and relies on
Python helper scripts found in the `bin/` directory.


This pipeline provides four main subcommands:

**Abundance:** Full workflow for read filtering, alignment, probability calculation, and abundance estimation.

**Build-database:** Build a custom GTDB-based database from input sequences and taxonomic mappings.

**Collapse-taxonomy:** Collapse a single-sample abundance table to a specified taxonomic rank.

**Combine-outputs:** Merge multiple relative-abundance tables into one multi-sample table at a given rank.

***Requirements***

- **Nextflow** 24.x or later
- **Docker**
- Optional environment variable `GTDB_DATABASE_DIR` pointing to a GTDB database
  (used as default for the `--db` parameter)




