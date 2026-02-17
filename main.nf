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
            // Directory mode (backward compatible)
            log.info "Discovering FASTQs from directory: ${params.input}"
            input_reads = channel
                .fromPath("${params.input}/**", checkIfExists: true)
                .filter { it.name =~ /\.(fastq|fq)(\.gz)?$/ }
                .ifEmpty { error "No FASTQ files found in ${params.input}" }
                .map { fastq -> tuple(fastq.simpleName.replaceAll(/\.(fastq|fq)$/, ''), fastq) }
        }

        ABUNDANCE(input_reads, fasta_file, tax_file)
    }

    // ─── SUBCOMMAND: BUILD-DATABASE ─────────────────────────────────────────
    else if (params.cmd == 'build-database') {
        if (!params.sequences || !params.seq2tax)
            error "For build-database you must provide --sequences and --seq2tax"

        BUILD_DATABASE(
            channel.value(file(params.sequences)),
            channel.value(file(params.seq2tax)),
            channel.value(file(params.taxonomy_list)),
            bin_ch
        )
    }

    // ─── SUBCOMMAND: COLLAPSE-TAXONOMY ──────────────────────────────────────
    else if (params.cmd == 'collapse-taxonomy') {
        if (!params.input_tsv)
            error "For collapse-taxonomy you must provide --input_tsv"

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





