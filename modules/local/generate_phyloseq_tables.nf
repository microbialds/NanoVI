// modules/local/generate_phyloseq_tables.nf — Consolidate per-sample outputs into phyloseq-compatible tables

process GENERATE_PHYLOSEQ_TABLES {

    label 'vi_output'

    input:
    path(sample_tsvs)

    output:
    path "otu_table_abundance.tsv", emit: otu_abundance
    path "otu_table_counts.tsv",    emit: otu_counts, optional: true   // only with --keep_counts
    path "tax_table.tsv",           emit: tax_table
    path "sample_metadata.tsv",     emit: sample_metadata
    path "versions.yml", emit: versions

    script:
    """
    generate_phyloseq_tables.py \
        ${sample_tsvs} \
        --output-dir .
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch otu_table_abundance.tsv otu_table_counts.tsv tax_table.tsv sample_metadata.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
