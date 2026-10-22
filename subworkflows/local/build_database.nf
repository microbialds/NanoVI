// subworkflows/local/build_database.nf — Build-database subworkflow

include { BUILD_DB } from '../../modules/local/build_database'

workflow BUILD_DATABASE {

    take:
    sequences       // path: input FASTA
    seq2tax         // path: seq2tax mapping
    taxonomy        // path: taxonomy list file or NCBI taxonomy directory
    taxonomy_flag   // val: '--taxonomy-list' or '--ncbi-taxonomy'

    main:
    BUILD_DB(sequences, seq2tax, taxonomy, taxonomy_flag)

    emit:
    versions = BUILD_DB.out.versions
}
