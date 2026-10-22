process GET_CIGAR {

    label 'vi_stats'

    input:
    tuple val(sample_id), path(sam_file)

    output:
  
    tuple val(sample_id), path("${sample_id}_cigar_info.json"), emit: cigar_json
    path "versions.yml", emit: versions

    script:
    """
    cigar_profile.py \\
        --sam "${sam_file}" \\
        --output "${sample_id}_cigar_info.json"
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    echo '{"log_probs": [], "zero_locs": [], "longest_align": {}}' > ${sample_id}_cigar_info.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}


