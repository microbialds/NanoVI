// modules/build_database/main.nf

process BUILD_DB {
    label 'vi_database'

    input:
    path sequences       // staged from sequences_ch
    path seq2tax         // staged from seq2tax_ch
    path taxonomy        // --taxonomy_list file or --ncbi_taxonomy directory
    val  taxonomy_flag   // '--taxonomy-list' or '--ncbi-taxonomy', matching taxonomy
    path copy_source, stageAs: 'copy_source/*'  // copy-number table or rrnDB pan-taxa file, or [] (F58); own directory, so its name cannot collide
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



// A database from a GTDB SSU file (F62, D44): bin/build_gtdb_db.py reads the lineages
// from the FASTA headers, so no sequence-to-taxonomy table or taxonomy is needed
process BUILD_GTDB_DB {
    label 'vi_database'

    input:
    path ssu                                         // --gtdb_ssu
    path metadata, stageAs: 'gtdb_metadata/*'        // --gtdb_metadata files, or [] (F58)
    path copy_table, stageAs: 'copy_source/*'        // --copy_numbers table, or []

    output:
    path "${params.db_name}/**"
    path "versions.yml", emit: versions

    script:
    def metadata_opt = metadata ? "--gtdb-metadata ${metadata}" : ''
    def copy_opt     = copy_table ? "--copy-numbers ${copy_table}" : ''
    """
    build_gtdb_db.py \\
        --ssu        ${ssu} \\
        --db-name    ${params.db_name} \\
        --output-dir ${params.db_name} \\
        --min-length ${params.gtdb_min_length} \\
        --cross-genus-rule ${params.cross_genus_rule} \\
        ${params.merge_gtdb_suffixes ? '--merge-gtdb-suffixes' : ''} \\
        ${metadata_opt} ${copy_opt}
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    mkdir -p ${params.db_name}
    touch ${params.db_name}/species_taxid.fasta ${params.db_name}/taxonomy.tsv
    touch ${params.db_name}/identical_16s.tsv ${params.db_name}/cross_genus_identical.tsv
    ${metadata || copy_table ? "touch ${params.db_name}/copy_numbers.tsv" : ''}
    ${params.merge_gtdb_suffixes ? "touch ${params.db_name}/merged_species.tsv" : ''}
    if [ "${params.cross_genus_rule}" != keep ]; then
        touch ${params.db_name}/cross_genus_removed.tsv ${params.db_name}/cross_genus_emptied.tsv
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
