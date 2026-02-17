process RUN_VI {

    label 'vi_inference'

    input:
    tuple val(sample_id), path(logp_json)

    output:
    tuple val(sample_id), path("${sample_id}_abundance.json"), emit: abundance

    script:
    """
    python3 -c '
import sys, json
sys.path.insert(0, "${projectDir}/bin")
from variational_inference import variational_inference_iterations

with open("${logp_json}") as f:
    data = json.load(f)

log_p_rgs = data["log_p_rgs"]
assigned_count = data["assigned_count"]
unassigned_count = data["unassigned_count"]

db_ids = list({tid for read in log_p_rgs.values() for tid in read[0]})

freq_full, freq_thresh, _ = variational_inference_iterations(
    log_p_rgs, db_ids, lli_thresh=0.0001, input_threshold=0.001
)

with open("${sample_id}_abundance.json", "w") as f:
    json.dump(freq_thresh if freq_thresh else freq_full, f)
'
    """
}

