process COLLAPSE {
    label 'vi_collapse'

    input:
      path input_path
      val  rank

    output:
      path "*-${rank}.tsv", emit: collapsed
      path "versions.yml", emit: versions

    script:
    """
    main.py collapse-taxonomy ${input_path} ${rank}
    python_versions.py "${task.process}" > versions.yml
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

