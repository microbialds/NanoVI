process INDEX_DB {

    label 'alignment'

    input:
    path species_fasta

    output:
    path "gtdb_index.mmi", emit: index

    script:
    """
    minimap2 -k ${params.kmer_size} -d gtdb_index.mmi ${species_fasta}
    """
}
