// subworkflows/local/build_database.nf — Build-database subworkflow

include { BUILD_DB } from '../../modules/local/build_database'

workflow BUILD_DATABASE {

    take:
    sequences       // path: input FASTA
    seq2tax         // path: seq2tax mapping
    taxonomy_list   // path: taxonomy list
    bin_dir         // path: bin/ directory

    main:
    BUILD_DB(sequences, seq2tax, taxonomy_list, bin_dir)
}
