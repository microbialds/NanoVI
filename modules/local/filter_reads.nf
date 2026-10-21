process FILTER_READS {

    label 'filter_reads'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: filtered

    // Additional outputs to publish
    path "fastplong/*.json.gz"
    path "fastplong/*.html.gz"
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p fastplong
    
    fastplong \\
        -i ${fastq_file} \\
        -o ${sample_id}_filtered.fastq \\
        --html fastplong/${sample_id}_report.html \\
        --json fastplong/${sample_id}_report.json \\
        --length_required ${params.min_length} \\
        --length_limit ${params.max_length} \\
        --disable_quality_filtering

    gzip -f ${sample_id}_filtered.fastq
    gzip -f fastplong/${sample_id}_report.html
    gzip -f fastplong/${sample_id}_report.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | sed 's/^fastplong //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p fastplong
    touch ${sample_id}_filtered.fastq.gz
    touch fastplong/${sample_id}_report.html.gz
    touch fastplong/${sample_id}_report.json.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: stub
    END_VERSIONS
    """
}

