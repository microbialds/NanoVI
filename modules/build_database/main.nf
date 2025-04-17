// modules/build_database/main.nf
nextflow.enable.dsl = 2

process BUILD_DB {
    label 'vi_database'

    // Publica todo lo que haya generado dentro de ${params.db_name}
    publishDir "${params.output_dir}/${params.db_name}", mode: 'copy'

    input:
    path sequences       // staged from sequences_ch
    path seq2tax         // staged from seq2tax_ch
    path taxonomy_list   // staged from taxonomy_list_ch

    script:
    """
    # 1) Copia tu carpeta bin/ al workdir
    cp -r ${projectDir}/bin bin

    # 2) Llama al subcomando con las rutas ya staged
    python3 bin/main.py build-database \\
        --sequences      ${sequences} \\
        --seq2tax        ${seq2tax} \\
        --taxonomy-list  ${taxonomy_list} \\
        ${params.db_name}
    """
}


