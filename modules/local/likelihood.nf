process LIKELIHOOD {

    label 'vi_stats'

    input:
    tuple val(sample_id), path(alignments)

    output:
    tuple val(sample_id), path("${sample_id}_read_scores.npz"), emit: scores
    tuple val(sample_id), path("${sample_id}_read_counts.json"), emit: counts
    path "versions.yml", emit: versions

    script:
    """
    read_likelihood.py \\
        --alignments "${alignments}" \\
        --scores "${sample_id}_read_scores.npz" \\
        --counts "${sample_id}_read_counts.json"
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch ${sample_id}_read_scores.npz
    echo '{"assigned_count": 0, "unassigned_count": 0}' > ${sample_id}_read_counts.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}
