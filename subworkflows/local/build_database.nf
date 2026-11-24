// subworkflows/local/build_database.nf — Build-database subworkflow

include { BUILD_DB; BUILD_GTDB_DB } from '../../modules/local/build_database'

workflow BUILD_DATABASE {

    take:
    sequences       // path: input FASTA
    seq2tax         // path: seq2tax mapping
    taxonomy        // path: taxonomy list file or NCBI taxonomy directory
    taxonomy_flag   // val: '--taxonomy-list' or '--ncbi-taxonomy'
    copy_source     // path: copy-number table or rrnDB pan-taxa file, or [] for none (F58)
    copy_flag       // val: '--copy-numbers', '--rrndb-pantaxa' or ''

    main:
    BUILD_DB(sequences, seq2tax, taxonomy, taxonomy_flag, copy_source, copy_flag)

    emit:
    versions = BUILD_DB.out.versions
}

workflow BUILD_GTDB_DATABASE {

    take:
    ssu             // path: GTDB SSU FASTA (F62)
    metadata        // path: GTDB bac120/ar53 metadata files, or [] for none (F58)
    copy_table      // path: copy-number table, or [] for none

    main:
    BUILD_GTDB_DB(ssu, metadata, copy_table)

    emit:
    versions = BUILD_GTDB_DB.out.versions
}
