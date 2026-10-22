process INDEX_DB {

    label 'alignment'

    // Built once per database, preset, k-mer size and minimap2 version (F26): index_name
    // records the first three and task.ext.minimap2_version (conf/base.config) the last,
    // so a later run finds the index in the cache and skips this step
    storeDir "${params.index_cache_dir ?: "${params.output_dir}/index_cache"}"

    input:
    path species_fasta
    val index_name

    output:
    path "${index_name}.minimap2-${task.ext.minimap2_version}.mmi", emit: index
    path "${index_name}.minimap2-${task.ext.minimap2_version}.versions.yml", emit: versions

    script:
    def index = "${index_name}.minimap2-${task.ext.minimap2_version}"
    """
    # k, w and homopolymer compression are fixed in the index, so it is built with
    # the alignment preset and k-mer size (F19); minimap2 applies -x before -k
    minimap2 -x ${params.type} -k ${params.kmer_size} -d ${index}.mmi ${species_fasta}
    cat <<-END_VERSIONS > ${index}.versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version)"
    END_VERSIONS
    """

    stub:
    def index = "${index_name}.minimap2-${task.ext.minimap2_version}"
    """
    touch ${index}.mmi
    cat <<-END_VERSIONS > ${index}.versions.yml
    "${task.process}":
        minimap2: stub
    END_VERSIONS
    """
}
