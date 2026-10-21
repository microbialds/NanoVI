process RUN_VI {

    label 'vi_inference'

    input:
    tuple val(sample_id), path(logp_json)

    output:
    tuple val(sample_id), path("${sample_id}_abundance.json"), emit: abundance
    path "versions.yml", emit: versions

    script:
    """
    estimate_abundance.py \\
        --logp "${logp_json}" \\
        --output "${sample_id}_abundance.json"
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    echo '{"abundance": {}, "below_detection_count": 0}' > ${sample_id}_abundance.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

