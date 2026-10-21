process ALIGN_READS {

    label 'alignment'

    input:
    tuple val(sample_id), path(filtered_fastq)
    path db_index

    output:
    tuple val(sample_id), path("${sample_id}_vi_alignments.sam"), emit: sam
    path "versions.yml", emit: versions

    script:
    """
    echo "Aligning ${sample_id} with minimap2..."

    minimap2 -ax ${params.type} \\
             -t ${task.cpus} \\
             -N ${params.N} \\
             -p 0.9 \\
             -K ${params.K} \\
             -k ${params.kmer_size} \\
             --split-prefix ${params.split_prefix} \\
             ${db_index} \\
             ${filtered_fastq} \\
             -o ${sample_id}_vi_alignments.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_vi_alignments.sam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
    END_VERSIONS
    """
}
