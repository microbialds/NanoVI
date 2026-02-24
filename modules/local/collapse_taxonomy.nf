process COLLAPSE {
    label 'vi_collapse'

    input:
      path input_path
      val  rank
      path bin_dir     // staged bin/ directory from bin_ch

    output:
      path "*-${rank}.tsv", emit: collapsed

    script:
    """
    python3 ${bin_dir}/main.py collapse-taxonomy ${input_path} ${rank}
    """
}

