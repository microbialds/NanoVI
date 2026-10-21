process GET_CIGAR {

    label 'vi_stats'

    input:
    tuple val(sample_id), path(sam_file)
    path bin_dir

    output:
  
    tuple val(sample_id), path("${sample_id}_cigar_info.json"), emit: cigar_json
    path "versions.yml", emit: versions

    script:
    """
    python3 -c '
import sys, json
sys.path.insert(0, "${bin_dir}")
from likelihood import get_cigar_op_log_probabilities

log_probs, zero_locs, longest_align = get_cigar_op_log_probabilities("${sam_file}", ${task.cpus})

with open("${sample_id}_cigar_info.json", "w") as f:
    json.dump({
        "log_probs": log_probs,
        "zero_locs": zero_locs,
        "longest_align": longest_align
    }, f)
'
    python3 ${bin_dir}/python_versions.py "${task.process}" > versions.yml
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


