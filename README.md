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
```

- [**Nextflow**](https://www.nextflow.io/docs/latest/install.html) v24.x or later
- [**Docker**](https://docs.docker.com/engine/install/)
- Optional environment variable `GTDB_DATABASE_DIR` pointing to a GTDB database
  (used as default for the `--db` parameter)



**2. Input file**

NanoVI accepts input files with the extensions .fastq or .fastq.gz. If you need to process multiple samples simultaneously, you should provide the path to the input folder.


**3. Pipeline usage**
```
nextflow run main.nf --cmd <subcommand> [parameters...]

```

Example commands:

Abundance estimation

```
nextflow run main.nf \
  --cmd abundance \
  --input /path/to/fastq \
  --output_dir /path/to/output_folder \
  --db /path/to/db \
  --taxonomy_tsv /path/to/db/taxonomy.tsv \
  --threads <number of threads>

```
Build database

```
nextflow run main.nf \
  --cmd build-database \
  --sequences /path/to/sequences.fasta \
  --seq2tax /path/to/seq2tax.tsv \
  --db_name custom_gtdb_database \
  --output_dir /path/to/output_db

```
Collapse taxonomy

```
nextflow run main.nf \
  --cmd collapse-taxonomy \
  --input_tsv /path/to/abundance_table.tsv \
  --rank genus

```
Combine outputs

```
nextflow run main.nf \
  --cmd combine-outputs \
  --input_dir /path/to/folder_with_abundance_tables \
  --rank species \
  --output_dir /path/to/combined_output

```

**4. Pipeline parameters (detailed)**
Parameter	Description	Default value
--cmd	Subcommand to run: abundance, build-database, collapse-taxonomy, or combine-outputs	'abundance'
--input	Path to input FASTQ file or directory	'data/'
--output_dir	Path to output directory	'results/'
--db	Path to GTDB database folder	'db/'
--taxonomy_tsv	Path to taxonomy TSV file	'db/taxonomy.tsv'
--threads	Number of CPU threads to use	8
--kmer_size	K-mer size for analysis	27
--N	Number of reads to sample	20
--K	Max number of kmers	1_000_000_000
--type	Sequencing type (e.g. 'map-ont' for Nanopore reads)	'map-ont'
--split_prefix	Prefix for temporary split files	'temp'
--min_length	Minimum read length to consider	500
--max_length	Maximum read length to consider	2000
--keep_counts	Whether to keep raw read counts (true or false)	false
--keep_files	Whether to keep intermediate files (true or false)


**4. Command example for Abundance Estimation**




