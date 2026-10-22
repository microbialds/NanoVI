// modules/combine_outputs/main.nf

process COMBINE {
    label 'vi_combine'

    input:
      path input_dir
      val  rank

    output:
      path "vi-combined-${rank}*.tsv", emit: combined
      path "versions.yml", emit: versions

    script:
    """
    mkdir -p input_data
    cp -v ${input_dir}/*_rel-abundance.tsv input_data/ || cp -v ${input_dir}/*_abundance.tsv input_data/
    main.py combine-outputs \\
        input_data \\
        ${rank}
    mv input_data/vi-combined-${rank}*.tsv ./
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch vi-combined-${rank}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}




