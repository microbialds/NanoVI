nextflow.enable.dsl = 2

process COLLAPSE {
    label 'vi_collapse'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
      path input_path
      val  rank

    output:
      path "*-${rank}.tsv", emit: collapsed

    script:
    """
    cp -r ${projectDir}/bin bin
    python3 bin/main.py collapse-taxonomy ${input_path} ${rank}
    """
}

