# NanoVI Output Documentation

## Abundance Subcommand

### Per-sample abundance table

`<sample>_rel-abundance.tsv`

A tab-separated file with columns:
- `lineage`: Full GTDB taxonomic lineage (domain;phylum;class;order;family;genus;species)
- `relative_abundance`: Estimated relative abundance (0.0 to 1.0)
- `count` (optional, with `--keep_counts`): Estimated number of reads assigned

### QC reports

`qc/<sample>_report.html.gz` and `qc/<sample>_report.json.gz`

FastpLong quality control reports showing read length distributions, quality scores, and filtering statistics.

### Intermediate files (with --keep_files)

- `filtered/<sample>_filtered.fastq.gz`: Length-filtered reads
- `sam/<sample>.sam`: Minimap2 alignments

## Build-Database Subcommand

Outputs the custom database directory with indexed reference sequences and taxonomy mappings.

## Collapse-Taxonomy Subcommand

`<sample>-<rank>.tsv`: Abundance table collapsed to the specified taxonomic rank.

## Combine-Outputs Subcommand

`vi-combined-<rank>.tsv`: Combined abundance matrix with samples as columns.

## Pipeline Information

- `pipeline_info/timeline.html`: Visual timeline of process execution
- `pipeline_info/report.html`: Nextflow execution statistics and resource usage
- `pipeline_info/trace.txt`: Tab-separated trace with per-task metrics
