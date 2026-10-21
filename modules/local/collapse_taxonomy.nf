process COLLAPSE {
    label 'vi_collapse'

    input:
      path input_path
      val  rank
      path bin_dir     // staged bin/ directory from bin_ch

    output:
      path "*-${rank}.tsv", emit: collapsed
      path "versions.yml", emit: versions

    script:
    """
    python3 ${bin_dir}/main.py collapse-taxonomy ${input_path} ${rank}
    python3 ${bin_dir}/python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch stub-${rank}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

