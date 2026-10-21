process WRITE_OUTPUT {

    label 'vi_output'

    input:
    tuple val(sample_id), path(abundance_json), path(logp_json)
    path taxonomy_tsv
    path bin_dir

    output:
    path "*.tsv", emit: result_tsv
    path "versions.yml", emit: versions

    script:
    """
    python3 -c '
import sys, json, pandas as pd
sys.path.insert(0, "${bin_dir}")
from taxonomy import freq_to_lineage_df

with open("${abundance_json}") as f:
    abundance_data = json.load(f)

freq = abundance_data["abundance"]
below_detection_count = abundance_data["below_detection_count"]

with open("${logp_json}") as f:
    logp_data = json.load(f)

assigned_count = logp_data["assigned_count"]
unassigned_count = logp_data["unassigned_count"]

taxonomy_df = pd.read_csv("${taxonomy_tsv}", sep="\\t", dtype=str).set_index("tax_id")

freq_to_lineage_df(
    freq,
    tsv_output_path="${sample_id}_rel-abundance",
    taxonomy_df=taxonomy_df,
    assigned_count=assigned_count,
    unassigned_count=unassigned_count,
    counts=${params.keep_counts ? 'True' : 'False'},
    below_detection_count=below_detection_count
)
    '
    python3 ${bin_dir}/python_versions.py "${task.process}" > versions.yml
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

