process FILTER_READS {

    label 'filter_reads'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered

    // Additional outputs to publish
    path "fastplong/*.json.gz"
    path "fastplong/*.html.gz"

    script:
    """
    mkdir -p fastplong
    
    fastplong \\
        -i ${fastq_file} \\
        -o ${sample_id}_filtered.fastq \\
        --html fastplong/${sample_id}_report.html \\
        --json fastplong/${sample_id}_report.json \\
        --length_required ${params.min_length} \\
        --length_limit ${params.max_length}

    gzip -f ${sample_id}_filtered.fastq
    gzip -f fastplong/${sample_id}_report.html
    gzip -f fastplong/${sample_id}_report.json
    """
}

