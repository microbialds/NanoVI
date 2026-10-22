process ALIGN_READS {

    label 'alignment'

    input:
    tuple val(sample_id), path(filtered_fastq)
    path db_index

    output:
    tuple val(sample_id), path("${sample_id}_vi_alignments.bam"), emit: bam
    path "versions.yml", emit: versions

    // minimap2's SAM goes straight into samtools as BAM (F27); the shell runs with
    // pipefail (conf/base.config), so a minimap2 failure fails the task. When minimap2
    // fails, the task exits with its status, not samtools': a memory kill (137) must
    // reach the retry in conf/base.config even if samtools then fails on the cut-off SAM
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
        | samtools view -b -o ${sample_id}_vi_alignments.bam - \\
        || { status=(\${PIPESTATUS[@]}); exit \$(( status[0] != 0 ? status[0] : status[1] )); }
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: "\$(minimap2 --version)"
        samtools: "\$(samtools --version | sed -n '1s/^samtools //p')"
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_vi_alignments.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
        samtools: stub
    END_VERSIONS
    """
}
