// modules/local/generate_phyloseq_tables.nf — Consolidate per-sample outputs into phyloseq-compatible tables

process GENERATE_PHYLOSEQ_TABLES {

    label 'vi_output'

    input:
    path(sample_tsvs)
    path bin_dir

    output:
    path "otu_table_abundance.tsv", emit: otu_abundance
    path "otu_table_counts.tsv",    emit: otu_counts
    path "tax_table.tsv",           emit: tax_table
    path "sample_metadata.tsv",     emit: sample_metadata

    script:
    """
    python3 ${bin_dir}/generate_phyloseq_tables.py \
        ${sample_tsvs} \
        --output-dir .
    """
}
