process WRITE_OUTPUT {

    label 'vi_output'

    input:
    tuple val(sample_id), path(abundance_json), path(logp_json)
    path taxonomy_tsv

    output:
    path "*.tsv", emit: result_tsv
    path "versions.yml", emit: versions

    script:
    """
    write_abundance_table.py \\
        --abundance "${abundance_json}" \\
        --logp "${logp_json}" \\
        --taxonomy "${taxonomy_tsv}" \\
        --output-prefix "${sample_id}_rel-abundance" \\
        ${params.keep_counts ? '--counts' : ''}
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch ${sample_id}_rel-abundance.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

