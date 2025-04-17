# **NanoVI**

Nextflow DSL2 implementation for variational inference of metagenomic abundances

This pipeline provides four main subcommands:

**Abundance:** Full workflow for read filtering, alignment, probability calculation, and abundance estimation.

**Build-database:** Build a custom GTDB-based database from input sequences and taxonomic mappings.

**Collapse-taxonomy:** Collapse a single-sample abundance table to a specified taxonomic rank.

**Combine-outputs:** Merge multiple relative-abundance tables into one multi-sample table at a given rank.
