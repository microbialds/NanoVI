process COMPUTE_LOGP {

    label 'vi_stats'

    input:
    tuple val(sample_id), path(sam_file), path(cigar_json)

    output:
    tuple val(sample_id), path("${sample_id}_logp_data.json"), emit: logp
    path "versions.yml", emit: versions

    script:
    """
    read_scores.py \\
        --sam "${sam_file}" \\
        --cigar-json "${cigar_json}" \\
        --output "${sample_id}_logp_data.json"
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    echo '{"log_p_rgs": {}, "assigned_count": 0, "unassigned_count": 0}' > ${sample_id}_logp_data.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

