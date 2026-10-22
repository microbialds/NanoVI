// main.nf – NanoVI: Taxonomic classification via variational inference
nextflow.enable.dsl = 2

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

// ─── PARAMETER TYPES ────────────────────────────────────────────────────────
// Defaults are in nextflow.config and descriptions in nextflow_schema.json. The
// types make Nextflow convert values given on the command line, which arrive as
// strings (on Nextflow 26.04, `--keep_counts false` was the truthy string "false"),
// and reject values of the wrong type (F40)
params {
    // General
    cmd: String
    output_dir: String

    // Abundance
    input: String?
    db: String?
    taxonomy_tsv: String?
    kmer_size: Integer
    N: Integer
    K: String
    type: String
    split_prefix: String
    min_length: Integer
    max_length: Integer
    min_quality: Integer?
    keep_counts: Boolean
    keep_files: Boolean

    // Abundance estimation (RUN_VI)
    estimator: String
    min_abundance: Float
    vi_tolerance: Float
    vi_max_iterations: Integer
    alpha_prior: Float
    vi_elbo_threshold: Float
    credible_intervals: Boolean

    // Build database
    db_name: String
    sequences: String?
    seq2tax: String?
    ncbi_taxonomy: String?
    taxonomy_list: String?

    // Post-processing
    input_tsv: String?
    input_dir: String?
    rank: String?
}

// ─── SUBWORKFLOW INCLUDES ───────────────────────────────────────────────────
include { ABUNDANCE }          from './subworkflows/local/abundance'
include { BUILD_DATABASE }     from './subworkflows/local/build_database'
include { COLLAPSE_TAXONOMY }  from './subworkflows/local/postprocessing'
include { COMBINE_OUTPUTS }    from './subworkflows/local/postprocessing'

// ─── MAIN WORKFLOW ──────────────────────────────────────────────────────────
workflow {

    // ─── PARAMETER VALIDATION ───────────────────────────────────────────────
    // Types, ranges, enums and paths from nextflow_schema.json (F40)
    validateParameters()

    // ─── VERSION LOGGING ────────────────────────────────────────────────────
    log.info """
    ================================================
      NanoVI Pipeline v${workflow.manifest.version}
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

        if (params.credible_intervals && params.estimator != 'vb')
            error "--credible_intervals needs --estimator vb: the EM estimator gives no posterior"

        // Input: samplesheet or directory
        if (params.input.endsWith('.csv')) {
            // Samplesheet mode
            log.info "Reading samplesheet: ${params.input}"
            // Rows checked against assets/schema_input.json: name pattern, unique
            // IDs, FASTQ extension, and the file exists (local, relative or remote)
            def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
            input_reads = channel.fromList(rows)
                .map { meta, fastq -> tuple(meta.id, fastq) }
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
        versions = ABUNDANCE.out.versions
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
            taxonomy_flag
        )
        versions = BUILD_DATABASE.out.versions
    }

    // ─── SUBCOMMAND: COLLAPSE-TAXONOMY ──────────────────────────────────────
    else if (params.cmd == 'collapse-taxonomy') {
        if (!params.input_tsv)
            error "For collapse-taxonomy you must provide --input_tsv"
        validateRank('collapse-taxonomy')

        COLLAPSE_TAXONOMY(
            channel.value(file(params.input_tsv)),
            channel.value(params.rank)
        )
        versions = COLLAPSE_TAXONOMY.out.versions
    }

    // ─── SUBCOMMAND: COMBINE-OUTPUTS ────────────────────────────────────────
    else if (params.cmd == 'combine-outputs') {
        if (!params.input_dir)
            error "For combine-outputs you must provide --input_dir"
        validateRank('combine-outputs')

        COMBINE_OUTPUTS(
            channel.value(file(params.input_dir)),
            channel.value(params.rank)
        )
        versions = COMBINE_OUTPUTS.out.versions
    }

    // ─── UNKNOWN COMMAND ────────────────────────────────────────────────────
    else {
        error "Unrecognized command: ${params.cmd}. Use: abundance, build-database, collapse-taxonomy, or combine-outputs"
    }

    // ─── SOFTWARE VERSIONS ──────────────────────────────────────────────────
    softwareVersionsToYAML(versions)
        .collectFile(storeDir: "${params.output_dir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)
}

// ─── HELPERS ────────────────────────────────────────────────────────────────

// Software versions (nf-core style): each process writes a versions.yml keyed by its
// full name; keep the process name only, drop the copies that per-sample tasks write,
// and add the pipeline and Nextflow versions
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file.text).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${workflow.manifest.version}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

def softwareVersionsToYAML(ch_versions) {
    return ch_versions
        .map { yaml_file -> processVersionsFromYAML(yaml_file) }
        .unique()
        .mix(channel.of(workflowVersionToYAML()))
}

// --rank must name a column of the per-sample tables; validateParameters() checks
// the value against the schema's enum, this checks that it is given
def validateRank(cmd) {
    def ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    if (!params.rank)
        error "For ${cmd} you must provide --rank (one of: ${ranks.join(', ')})"
}
