process RUN_VI {

    label 'vi_inference'

    input:
    tuple val(sample_id), path(read_scores)

    output:
    tuple val(sample_id), path("${sample_id}_abundance.json"), emit: abundance
    path "versions.yml", emit: versions

    script:
    """
    estimate_abundance.py \\
        --scores "${read_scores}" \\
        --output "${sample_id}_abundance.json" \\
        --estimator ${params.estimator} \\
        --min-abundance ${params.min_abundance} \\
        --tolerance ${params.vi_tolerance} \\
        --max-iterations ${params.vi_max_iterations} \\
        --alpha-prior ${params.alpha_prior} \\
        --elbo-threshold ${params.vi_elbo_threshold}
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    echo '{"abundance": {}, "below_detection_count": 0, "read_summary": {}}' > ${sample_id}_abundance.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

