process INDEX_DB {

    label 'alignment'

    input:
    path species_fasta

    output:
    path "gtdb_index.mmi", emit: index
    path "versions.yml", emit: versions

    script:
    """
    # k, w and homopolymer compression are fixed in the index, so it is built with
    # the alignment preset and k-mer size (F19); minimap2 applies -x before -k
    minimap2 -x ${params.type} -k ${params.kmer_size} -d gtdb_index.mmi ${species_fasta}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version)"
    END_VERSIONS
    """

    stub:
    """
    touch gtdb_index.mmi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
    END_VERSIONS
    """
}
