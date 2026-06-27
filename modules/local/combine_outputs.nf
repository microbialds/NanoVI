// modules/combine_outputs/main.nf

process COMBINE {
    label 'vi_combine'

    input:
      path input_dir
      val  rank
      path bin_dir     // staged bin/ directory from bin_ch

    output:
      path "vi-combined-${rank}*.tsv", emit: combined

    script:
    """
    mkdir -p input_data
    cp -v ${input_dir}/*_rel-abundance.tsv input_data/ || cp -v ${input_dir}/*_abundance.tsv input_data/
    python3 ${bin_dir}/main.py combine-outputs \\
        input_data \\
        ${rank}
    mv input_data/vi-combined-${rank}*.tsv ./
    """

    stub:
    """
    touch vi-combined-${rank}.tsv
    """
}




