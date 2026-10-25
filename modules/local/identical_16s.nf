process IDENTICAL_16S {

    label 'vi_database'

    // Built once per database FASTA, taxonomy and identical_16s.py (D24), when the database directory
    // has no identical_16s.tsv; kept beside the cached indexes, so a later run finds
    // it and skips this step. Every output goes into the shared cache, so the
    // versions file carries the table's name
    storeDir "${params.index_cache_dir ?: "${params.output_dir}/index_cache"}"

    input:
    path species_fasta
    path taxonomy_tsv
    val table_name

    output:
    path "${table_name}.tsv", emit: table
    path "${table_name}.versions.yml", emit: versions

    script:
    """
    identical_16s.py \\
        --fasta "${species_fasta}" \\
        --taxonomy "${taxonomy_tsv}" \\
        --output "${table_name}.tsv"
    python_versions.py "${task.process}" > ${table_name}.versions.yml
    """

    stub:
    """
    printf 'tax_id\\tspecies\\tgenus\\tgroup\\tgroup_size\\tidentical_16s_species\\n' > ${table_name}.tsv
    cat <<-END_VERSIONS > ${table_name}.versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
