// modules/combine_outputs/main.nf
nextflow.enable.dsl = 2

process COMBINE {
    label 'vi_combine'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
      path input_dir
      val  rank

    output:
      path "vi-combined-${rank}*.tsv", emit: combined

    script:
    """
    cp -r ${projectDir}/bin bin
    mkdir -p input_data
    cp -v ${input_dir}/*_rel-abundance.tsv input_data/ || cp -v ${input_dir}/*_abundance.tsv input_data/
    python3 bin/main.py combine-outputs \\
        input_data \\
        ${rank}
    mv input_data/vi-combined-${rank}*.tsv ./
    """
}




