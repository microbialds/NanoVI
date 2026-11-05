// modules/build_database/main.nf

process BUILD_DB {
    label 'vi_database'

    input:
    path sequences       // staged from sequences_ch
    path seq2tax         // staged from seq2tax_ch
    path taxonomy        // --taxonomy_list file or --ncbi_taxonomy directory
    val  taxonomy_flag   // '--taxonomy-list' or '--ncbi-taxonomy', matching taxonomy
    path copy_source     // copy-number table or rrnDB pan-taxa file, or [] (F58)
    val  copy_flag       // '--copy-numbers', '--rrndb-pantaxa' or ''

    output:
    path "${params.db_name}/**"
    path "versions.yml", emit: versions

    script:
    """
    main.py build-database \\
        --sequences      ${sequences} \\
        --seq2tax        ${seq2tax} \\
        ${taxonomy_flag} ${taxonomy} \\
        --cross-genus-rule ${params.cross_genus_rule} \\
        ${copy_flag ? "${copy_flag} ${copy_source}" : ''} \\
        ${params.db_name}
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    mkdir -p ${params.db_name}
    touch ${params.db_name}/species_taxid.fasta ${params.db_name}/taxonomy.tsv
    touch ${params.db_name}/identical_16s.tsv ${params.db_name}/cross_genus_identical.tsv
    ${copy_flag ? "touch ${params.db_name}/copy_numbers.tsv" : ''}
    if [ "${params.cross_genus_rule}" != keep ]; then
        touch ${params.db_name}/cross_genus_removed.tsv ${params.db_name}/cross_genus_emptied.tsv
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}


