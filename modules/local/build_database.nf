// modules/build_database/main.nf

process BUILD_DB {
    label 'vi_database'

    input:
    path sequences       // staged from sequences_ch
    path seq2tax         // staged from seq2tax_ch
    path taxonomy_list   // staged from taxonomy_list_ch
    path bin_dir         // staged bin/ directory from bin_ch

    output:
    path "${params.db_name}/**"

    script:
    """
    python3 ${bin_dir}/main.py build-database \\
        --sequences      ${sequences} \\
        --seq2tax        ${seq2tax} \\
        --taxonomy-list  ${taxonomy_list} \\
        ${params.db_name}
    """
}


