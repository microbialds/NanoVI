process WRITE_OUTPUT {

    label 'vi_output'

    input:
    tuple val(sample_id), path(abundance_json), path(read_counts)
    path taxonomy_tsv
    path identical_16s   // the database's identical_16s.tsv, or [] when there is none (D24)

    output:
    path "${sample_id}_rel-abundance.tsv", emit: result_tsv
    path "${sample_id}_groups.tsv", emit: groups, optional: !params.ambiguity_groups
    path "versions.yml", emit: versions

    script:
    // Without the table the identical_16s_species column is written empty
    def identical = identical_16s ? "\"${identical_16s}\"" : "''"
    """
    write_abundance_table.py \\
        --abundance "${abundance_json}" \\
        --read-counts "${read_counts}" \\
        --taxonomy "${taxonomy_tsv}" \\
        --output-prefix "${sample_id}_rel-abundance" \\
        --identical-16s ${identical} \\
        ${params.binomial_names ? '--binomial-names' : ''} \\
        ${params.credible_intervals ? '--credible-intervals' : ''} \\
        ${params.keep_counts ? '--counts' : ''} \\
        ${params.ambiguity_groups ? "--groups-output ${sample_id}_groups" : ''}
    python_versions.py "${task.process}" > versions.yml
    """

    stub:
    """
    touch ${sample_id}_rel-abundance.tsv
    ${params.ambiguity_groups ? "touch ${sample_id}_groups.tsv" : ''}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    """
}

