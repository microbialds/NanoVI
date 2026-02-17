process ALIGN_READS {

    label 'alignment'

    input:
    tuple val(sample_id), path(filtered_fastq)
    path db_index

    output:
    tuple val(sample_id), path("${sample_id}_vi_alignments.sam"), emit: sam

    publishDir "${params.output_dir}/sam", mode: 'copy', enabled: params.keep_files

    script:
    """
    echo "🧬 Alineando ${sample_id} con minimap2..."

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
    """
}
