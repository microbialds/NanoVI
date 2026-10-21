// main.nf – NanoVI: Taxonomic classification via variational inference
nextflow.enable.dsl = 2

// ─── SUBWORKFLOW INCLUDES ───────────────────────────────────────────────────
include { ABUNDANCE }          from './subworkflows/local/abundance'
include { BUILD_DATABASE }     from './subworkflows/local/build_database'
include { COLLAPSE_TAXONOMY }  from './subworkflows/local/postprocessing'
include { COMBINE_OUTPUTS }    from './subworkflows/local/postprocessing'

// ─── MAIN WORKFLOW ──────────────────────────────────────────────────────────
workflow {

    // ─── BIN DIRECTORY CHANNEL ──────────────────────────────────────────────
    bin_ch = channel.value(file("${projectDir}/bin"))

    // ─── VERSION LOGGING ────────────────────────────────────────────────────
    log.info """
    ================================================
      NanoVI Pipeline v1.0
    ================================================
      Command    : ${params.cmd}
      Input      : ${params.input}
      Database   : ${params.db}
      Output     : ${params.output_dir}
      K-mer size : ${params.kmer_size}
    ================================================
    """.stripIndent()

    // ─── SUBCOMMAND: ABUNDANCE ──────────────────────────────────────────────
    if (params.cmd == 'abundance') {

        // Validate inputs
        if (!params.input)
            error "Please provide --input (samplesheet CSV or FASTQ directory)"
        if (!params.db)
            error "Please provide --db (path to reference database)"

        def db_dir = file(params.db)
        if (!db_dir.exists())
            error "Database directory not found: ${params.db}"

        def fasta_file = file("${params.db}/species_taxid.fasta")
        if (!fasta_file.exists())
            error "Reference FASTA not found: ${params.db}/species_taxid.fasta"

        def taxonomy_tsv = params.taxonomy_tsv ?: "${params.db}/taxonomy.tsv"
        def tax_file = file(taxonomy_tsv)
        if (!tax_file.exists())
            error "Taxonomy file not found: ${taxonomy_tsv}"

        if (params.min_length >= params.max_length)
            error "min_length (${params.min_length}) must be less than max_length (${params.max_length})"

        // Input: samplesheet or directory
        if (params.input.endsWith('.csv')) {
            // Samplesheet mode
            log.info "Reading samplesheet: ${params.input}"
            def rows = new SamplesheetParser().parseSamplesheet(params.input)
            input_reads = channel.from(rows)
                .map { sample, fastq -> tuple(sample, file(fastq)) }
        } else {
            // Directory mode (backward compatible): the sample ID is the file name
            // without its FASTQ extension, so S1.run1.fastq.gz gives S1.run1
            log.info "Discovering FASTQs from directory: ${params.input}"
            def fastq_ext = /\.(fastq|fq)(\.gz)?$/
            def fastqs = files("${params.input}/**").findAll { f -> f.name =~ fastq_ext }
            if (!fastqs)
                error "No FASTQ files found in ${params.input}"
            def by_sample = fastqs.groupBy { f -> f.name.replaceAll(fastq_ext, '') }
            def clashes = by_sample.findAll { id, fs -> fs.size() > 1 }
            if (clashes) {
                def listing = clashes.collect { id, fs ->
                    "  ${id}: ${fs.collect { it.toString() }.sort().join(', ')}"
                }.sort().join('\n')
                error "Several FASTQ files in ${params.input} give the same sample ID; rename them or use a samplesheet:\n${listing}"
            }
            def bad_ids = by_sample.keySet().findAll { id -> !(id ==~ /[A-Za-z0-9._-]+/) }.sort()
            if (bad_ids)
                error "Sample IDs from FASTQ file names may only use letters, digits, '.', '_' and '-'; rename these files or use a samplesheet:\n" +
                    bad_ids.collect { id -> "  ${by_sample[id][0]}" }.join('\n')
            input_reads = channel.fromList(
                by_sample.collect { id, fs -> tuple(id, fs[0]) }.sort { it[0] }
            )
        }

        ABUNDANCE(input_reads, fasta_file, tax_file)
    }

    // ─── SUBCOMMAND: BUILD-DATABASE ─────────────────────────────────────────
    else if (params.cmd == 'build-database') {
        if (!params.sequences || !params.seq2tax)
            error "For build-database you must provide --sequences and --seq2tax"
        if (!params.taxonomy_list == !params.ncbi_taxonomy)
            error "For build-database provide exactly one of --taxonomy_list or --ncbi_taxonomy"

        // main.py takes the taxonomy source as --taxonomy-list <file> or --ncbi-taxonomy <dir>
        def taxonomy_flag = params.taxonomy_list ? '--taxonomy-list' : '--ncbi-taxonomy'
        def taxonomy_path = params.taxonomy_list ?: params.ncbi_taxonomy

        BUILD_DATABASE(
            channel.value(file(params.sequences)),
            channel.value(file(params.seq2tax)),
            channel.value(file(taxonomy_path)),
            taxonomy_flag,
            bin_ch
        )
    }

    // ─── SUBCOMMAND: COLLAPSE-TAXONOMY ──────────────────────────────────────
    else if (params.cmd == 'collapse-taxonomy') {
        if (!params.input_tsv)
            error "For collapse-taxonomy you must provide --input_tsv"
        validateRank('collapse-taxonomy')

        COLLAPSE_TAXONOMY(
            channel.value(file(params.input_tsv)),
            channel.value(params.rank),
            bin_ch
        )
    }

    // ─── SUBCOMMAND: COMBINE-OUTPUTS ────────────────────────────────────────
    else if (params.cmd == 'combine-outputs') {
        if (!params.input_dir)
            error "For combine-outputs you must provide --input_dir"
        validateRank('combine-outputs')

        COMBINE_OUTPUTS(
            channel.value(file(params.input_dir)),
            channel.value(params.rank),
            bin_ch
        )
    }

    // ─── UNKNOWN COMMAND ────────────────────────────────────────────────────
    else {
        error "Unrecognized command: ${params.cmd}. Use: abundance, build-database, collapse-taxonomy, or combine-outputs"
    }
}

// ─── HELPERS ────────────────────────────────────────────────────────────────

// --rank must name a column of the per-sample tables
def validateRank(cmd) {
    def ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    if (!params.rank)
        error "For ${cmd} you must provide --rank (one of: ${ranks.join(', ')})"
    if (!(params.rank in ranks))
        error "Invalid --rank '${params.rank}' for ${cmd}. Use one of: ${ranks.join(', ')}"
}
