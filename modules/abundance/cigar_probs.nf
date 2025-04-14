// modules/abundance/cigar_probs.nf

process GET_CIGAR {

    label 'vi_stats'

    input:
    path sam_file
    val threads

    output:
    path "cigar_info.json", emit: cigar_json

    script:
    """
    python3 -c '
import sys, json
sys.path.insert(0, "./bin")
from variational_inference import get_cigar_op_log_probabilities

log_probs, zero_locs, longest_align = get_cigar_op_log_probabilities("${sam_file}", ${threads})

with open("cigar_info.json", "w") as f:
    json.dump({
        "log_probs": log_probs,
        "zero_locs": zero_locs,
        "longest_align": longest_align
    }, f)
    '
    """
}

