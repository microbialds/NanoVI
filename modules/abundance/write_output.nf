// modules/abundance/write_output.nf

process WRITE_OUTPUT {

    label 'vi_output'

    input:
    path abundance_json
    path taxonomy_tsv
    val assigned_count
    val unassigned_count
    val output_prefix

    output:
    path "abundance_output.tsv", emit: result_tsv

    publishDir "${params.output_dir}", mode: 'copy'

    script:
    """
    python3 -c '
import sys, json, pandas as pd
sys.path.insert(0, "./bin")
from taxonomy import freq_to_lineage_df

# Cargar datos
with open("${abundance_json}") as f:
    freq = json.load(f)

taxonomy_df = pd.read_csv("${taxonomy_tsv}", sep="\\t", dtype=str).set_index("tax_id")

# Escribir resultado localmente
freq_to_lineage_df(
    freq,
    tsv_output_path="abundance_output",  # ⚠️ Sin rutas absolutas
    taxonomy_df=taxonomy_df,
    assigned_count=${assigned_count},
    unassigned_count=${unassigned_count},
    counts=${params.keep_counts ? 'True' : 'False'}

)
    '
    """
}

