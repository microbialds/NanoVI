// modules/build_database/main.nf

process BUILD_DB {
    label 'vi_database'

    input:
    path sequences       // staged from sequences_ch
    path seq2tax         // staged from seq2tax_ch
    path taxonomy        // --taxonomy_list file or --ncbi_taxonomy directory
    val  taxonomy_flag   // '--taxonomy-list' or '--ncbi-taxonomy', matching taxonomy
    path bin_dir         // staged bin/ directory from bin_ch

    output:
    path "${params.db_name}/**"
    path "versions.yml", emit: versions

    script:
    """
    python3 ${bin_dir}/main.py build-database \\
        --sequences      ${sequences} \\
        --seq2tax        ${seq2tax} \\
        ${taxonomy_flag} ${taxonomy} \\
        ${params.db_name}
    python3 ${bin_dir}/python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    mkdir -p ${params.db_name}
    touch ${params.db_name}/species_taxid.fasta ${params.db_name}/taxonomy.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}


