// modules/abundance/filter_reads.nf

process FILTER_READS {

    label 'filter_reads'

    input:
    path fastq_file

    output:
    tuple val("${fastq_file.simpleName}"), path("*.fastq.gz"), emit: filtered
    path "fastplong/*.json.gz"
    path "fastplong/*.html.gz"
    publishDir "${params.output_dir}", mode: 'copy'
    script:
    """
    sample_name=\$(basename ${fastq_file} | sed 's/.fastq.gz//; s/.fq.gz//; s/.fastq//; s/.fq//')

    mkdir -p fastplong

    fastplong \\
        -i ${fastq_file} \\
        -o \${sample_name}_filtered.fastq \\
        --html fastplong/\${sample_name}_report.html \\
        --json fastplong/\${sample_name}_report.json \\
        --length_required ${params.min_length} \\
        --length_limit ${params.max_length}

    gzip -f \${sample_name}_filtered.fastq
    gzip -f fastplong/\${sample_name}_report.html
    gzip -f fastplong/\${sample_name}_report.json
    """
}
