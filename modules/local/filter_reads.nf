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

    // fastplong writes the filtered reads gzipped (F27) with the CPUs the task asks for
    // (F28; its default is 3 worker threads); the output does not depend on the count
    script:
    // --min_quality turns on fastplong's quality filter on mean read quality (D5; off by default):
    // the per-base rule is switched off, fastplong's limit of 5 N bases per read still applies
    def quality_filter = params.min_quality != null ?
        "--mean_qual ${params.min_quality} --unqualified_percent_limit 100" :
        '--disable_quality_filtering'
    """
    mkdir -p fastplong

    fastplong \\
        -i ${fastq_file} \\
        -o ${sample_id}_filtered.fastq.gz \\
        --thread ${task.cpus} \\
        --html fastplong/${sample_id}_report.html \\
        --json fastplong/${sample_id}_report.json \\
        --length_required ${params.min_length} \\
        --length_limit ${params.max_length} \\
        ${quality_filter}

    gzip -f fastplong/${sample_id}_report.html
    gzip -f fastplong/${sample_id}_report.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: "\$(fastplong --version 2>&1 | sed 's/^fastplong //')"
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

