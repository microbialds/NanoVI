process COMPUTE_LOGP {

    label 'vi_stats'

    input:
    tuple val(sample_id), path(sam_file), path(cigar_json)

    output:
    tuple val(sample_id), path("${sample_id}_logp_data.json"), emit: logp

    script:
    """
    python3 -c '
import sys, json
sys.path.insert(0, "${projectDir}/bin")
from variational_inference import log_prob_rgs_dict

with open("${cigar_json}") as f:
    data = json.load(f)

log_probs     = data["log_probs"]
zero_locs     = data["zero_locs"]
longest_align = data["longest_align"]

log_p_rgs, unassigned_count, assigned_count = log_prob_rgs_dict(
    "${sam_file}", log_probs, longest_align, zero_locs
)

with open("${sample_id}_logp_data.json", "w") as f:
    json.dump({
        "log_p_rgs": log_p_rgs,
        "assigned_count": assigned_count,
        "unassigned_count": unassigned_count
    }, f)
'
    """
}

