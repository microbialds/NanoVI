# **NanoVI**

Pipeline for estimating relative abundances of Oxford Nanopore 16S reads using a variational inference approach. The workflow is written in Nextflow and relies on
Python helper scripts found in the `bin/` directory.


### 1. Pipeline structure

The NanoVI pipeline follows a modular architecture implemented in **Nextflow**.

- The main script `main.nf` acts as an entry point and dispatches the execution based on the `--cmd` parameter.
- Each subcommand (`abundance`, `build-database`, `collapse-taxonomy`, `combine-outputs`) is implemented as an independent workflow in the `modules/` directory.
- The `abundance` workflow is further divided into modular Nextflow processes (`filter_reads`, `alignment`, `inference`, etc.), allowing clear separation of steps and easier maintenance.
- Helper Python scripts used during processing are located in the `bin/` directory.
- A `test/` directory is provided with example data to test the pipeline.

This structure ensures that the pipeline is:

- Easy to extend with new subcommands or modules.
- Highly reproducible and portable through the use of **Nextflow** and **Docker**.
- User-friendly: each subcommand can be run independently by specifying the appropriate parameters.

Below is the current directory structure:

```
NanoVI/
├── main.nf                       # Main Nextflow script: dispatches to subcommand workflows
├── config/
│   ├── containers.config         # Docker container configuration
├── modules/                      # Modular Nextflow workflows by subcommand
│   ├── abundance/
│   │   ├── filter_reads.nf
│   │   ├── alignment.nf
│   │   ├── cigar_probs.nf
│   │   ├── log_prob_rgs.nf
│   │   ├── inference.nf
│   │   └── write_output.nf
│   ├── build_database/
│   │   └── main.nf
│   ├── collapse_taxonomy/
│   │   └── main.nf
│   └── combine_outputs/
│       └── main.nf
├── bin/                          # Python helper scripts
├── test/                      # Example FASTQ data for testing the pipeline
├── README.md                     # Project documentation
└── nextflow.config               # Pipeline configuration and default parameters

```


### 2. Installation

Clone this repository, and install:

```
git clone https://github.com/microbialds/NanoVI
cd NanoVI
```

- [**Nextflow**](https://www.nextflow.io/docs/latest/install.html) v24.x or later
- [**Docker**](https://docs.docker.com/engine/install/)
- Optional environment variable `GTDB_DATABASE_DIR` pointing to a GTDB database
  (used as default for the `--db` parameter)



### 3. Input file

NanoVI accepts input files with the extensions .fastq or .fastq.gz. If you need to process multiple samples simultaneously, you should provide the path to the input folder.


### 4. Pipeline usage

This pipeline provides four main subcommands:

 - Abundance: Full workflow for read filtering, alignment, probability calculation, and abundance estimation.

- Build-database: Build a custom GTDB-based database from input sequences and taxonomic mappings.

- Collapse-taxonomy: Collapse a single-sample abundance table to a specified taxonomic rank.

- Combine-outputs: Merge multiple relative-abundance tables into one multi-sample table at a given rank.


```
nextflow run main.nf --cmd <subcommand> [parameters...]

```


**Example commands:**

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


### 5. Pipeline parameters


- **Abundance parameters**

| Parameter                      | Default                     | Description                                                                                                                                                   |
|-------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--cmd`                        | `user must provide`        | Subcommand to run: `abundance`, `build-database`, `collapse-taxonomy`, or `combine-outputs`                                                                     |
| `--input`                      | `user must provide`        | Path to input FASTQ file or directory                                                                                                                          |
| `--output_dir`                  | `user must provide`       | Directory for output results                                                                                                                                  |
| `--db`                         | `$GTDB_DATABASE_DIR`        | Path to GTDB database folder                                                                                                                                  |
| `--taxonomy_tsv`               | `db/taxonomy.tsv`           | Path to taxonomy TSV file corresponding to database                                                                                                           |
| `--threads`                    | `8`                         | Number of CPU threads to use                                                                                                                                  |
| `--kmer_size`                  | `27`                        | K-mer size for minimap2                                                                                                                                       |
| `--N`                          | `20`                        | max number of alignments utilized for each read in minimap2                                                                                                                                 |
| `--K`                          | `1000000000`                | minibatch size for mapping in minimap2                                                                                                                         |
| `--type`                       | `map-ont`                   | Sequencing type (e.g. `'map-ont'` for Nanopore reads)                                                                                                          |
| `--split_prefix`               | `temp`                      | Prefix for temporary split files                                                                                                                              |
| `--min_length`                 | `500`                       | Minimum read length to fastplong filter                                                                                                                       |
| `--max_length`                 | `2000`                      | Maximum read length to fastplong filter                                                                                                                       |
| `--keep_counts`                | `FALSE`                     | Include estimated read counts for each species in output                                                                                                      |
| `--keep_files`                 | `FALSE`                     | Keep intermediate files in output directory (alignments [.sam], reads of specified length [.fa])                                                               |



- **Build database parameters**

| Parameter                      | Default                     | Description                                                                                                                                                   |
|-------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--db_name`                    | `db_custom`                 | Name of the custom database to build                                                                                                                           |
| `--sequences`                  | `null`                      | Input FASTA sequences in fasta format                                                                                                                                        |
| `--seq2tax`                    | `null`                      | Sequence-to-taxonomy mapping file in tsv format (seq2tax.map.tsv)                                                                                                                            |
| `--ncbi_taxonomy`              | `null`                      | NCBI taxonomy dump directory (names.dmp & nodes.dmp files)                                                                                                                                  |
| `--taxonomy_list`              | `null`                      | List of taxonomy terms to include in tsv format                                                                                                                            |



- **Collapse taxonomy parameters**

| Parameter                      | Default                     | Description                                                                                                                                                   |
|-------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input_tsv`                  | `null`                      | Input abundance table TSV from output directory                                                                                                                                     |
| `--rank`                       | `null`                      | Taxonomic rank to collapse to (e.g., `phylum`, `genus`, `species`)                                                                                              |



- **Combine outputs parameters**

| Parameter                      | Default                     | Description                                                                                                                                                   |
|-------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input_dir`                  | `null`                      | Directory containing abundance tables from output directory                                                                                                                        |
| `--rank`                       | `null`                      | Taxonomic rank to combine tables                                                                                                                           |


### 4. GTDB database

GTDB database version 220 available, in the following link:

https://drive.google.com/drive/folders/1VurYfWfksEuokUCQCCSEbqJ8jBAEUXAw?usp=sharing

You can also create your own database or use the EMU database provided by the EMU authors

### 5. Pipeline testing with example dataset

We include a dataset corresponding to a mock community for running the pipeline. This dataset is located in the test/ directory and serves as a benchmark for validating the pipeline's performance.







