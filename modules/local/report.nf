// modules/local/report.nf — The run report: an offline HTML site, a PDF and paper-ready figures (--report; F60, D42)

process REPORT {

    label 'report'

    input:
    // What the run publishes, staged at its published paths so that nanovi_report.py reads
    // a results directory (the same script renders one after the run)
    path(sample_tables, stageAs: 'results/sample_outputs/*')       // <sample>_rel-abundance.tsv, and the optional tables
    path(consolidated, stageAs: 'results/consolidated_output/*')   // otu_table_abundance, tax_table, sample_metadata
    path(fastplong_json, stageAs: 'results/qc/fastplong/*')        // <sample>_report.json.gz
    path(software_versions, stageAs: 'results/pipeline_info/*')    // software_versions.yml
    path(metadata)                                                 // --report_metadata, or []
    // assets/report and a hash of its files' contents: a directory input is cached by its
    // path, size and mtime, so without the hash an edited template would not rerun on -resume
    tuple path(report_assets, stageAs: 'report_assets'), val(report_assets_hash)
    val(run_params)   // the run's parameters as JSON, listed in the methods
    val(report_date)  // YYYY-MM-DD: the day the run started (UTC), so two renders give the same files

    output:
    path "report", emit: report
    path "versions.yml", emit: versions

    script:
    def metadata_arg = metadata ? "--metadata '${metadata}'" : ''
    """
    # matplotlib's and fontconfig's caches in the task directory (a read-only image, maybe no home)
    export MPLCONFIGDIR="\$PWD/.mplconfig" XDG_CACHE_HOME="\$PWD/.cache"
    cat > run_parameters.json <<'END_PARAMS'
    ${run_params}
    END_PARAMS
    nanovi_report.py \\
        --results results \\
        --outdir report \\
        --assets report_assets \\
        --params run_parameters.json \\
        --date ${report_date} \\
        ${metadata_arg}
    nanovi_report.py --versions "${task.process}" > versions.yml
    cp versions.yml report/versions.yml
    """

    stub:
    """
    mkdir -p report/figures report/assets
    touch report/index.html report/methods.html report/nanovi_report.pdf report/assets/site.css
    touch report/figures/read_fate.svg report/figures/read_fate.pdf report/figures/read_fate.png report/figures/read_fate.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: stub
    END_VERSIONS
    cp versions.yml report/versions.yml
    """
}
