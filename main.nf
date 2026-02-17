// main.nf – Nextflow version of the VI abundance pipeline
nextflow.enable.dsl = 2

// ─── INCLUDES (DSL2) ────────────────────────────────────────────────────────
include { FILTER_READS }  from './modules/local/filter_reads'
include { INDEX_DB }      from './modules/local/index_db'
include { ALIGN_READS }   from './modules/local/align_reads'
include { GET_CIGAR }     from './modules/local/cigar_probs'
include { COMPUTE_LOGP }  from './modules/local/log_prob_rgs'
include { RUN_VI }        from './modules/local/inference'
include { WRITE_OUTPUT }  from './modules/local/write_output'

include { BUILD_DB }      from './modules/local/build_database'
include { COLLAPSE }      from './modules/local/collapse_taxonomy'
include { COMBINE }       from './modules/local/combine_outputs'

workflow {

  // ─── SUBCOMAND: ABUNDANCE ────────────────────────────────────────────────
  if( params.cmd == 'abundance' ) {

    // Read FASTQ files
    input_reads = Channel
      .fromPath("${params.input}/**", checkIfExists: true)
      .filter { it.name.endsWith('.fastq') || it.name.endsWith('.fq') || it.name.endsWith('.fastq.gz') || it.name.endsWith('.fq.gz') }
      .ifEmpty { error "No FASTQ files found in ${params.input}" }

    // Step 1: Filter reads (tuple: sample_id, fastq)
    named_reads = FILTER_READS(input_reads).filtered

    // Step 2a: Build minimap2 index (once)
    fasta_file = file("${params.db}/species_taxid.fasta")
    db_index = INDEX_DB(fasta_file)

    // Step 2b: Alignment (uses pre-built index)
    sam_file = ALIGN_READS(
      named_reads,    // Tuple (sample_id, fastq.gz)
      db_index        // Pre-built minimap2 index
    )

    // Step 3: Get CIGAR probabilities
    cigar_info = GET_CIGAR(sam_file)

    // ✅ Join SAM and CIGAR for sample_id
    sam_and_cigar = sam_file.join(cigar_info)

    // Step 4: Log-probabilidades calculate
    logp_data = COMPUTE_LOGP(sam_and_cigar)


    // Step 5: Variational inference algorithm
    freq_output = RUN_VI(logp_data)

    // Step 6+7: Join abundance with logp data and write output
    vi_with_counts = freq_output.join(logp_data)
    // vi_with_counts is: tuple(sample_id, abundance.json, logp_data.json)

    WRITE_OUTPUT(
      vi_with_counts,
      file(params.taxonomy_tsv)
    )

   }

  // ─── SUBCOMAND: BUILD-DATABASE ───────────────────────────────────────────
  else if( params.cmd == 'build-database' ) {
    if( !params.sequences || !params.seq2tax )
      error "For build-database you must give --sequences and --seq2tax"

    BUILD_DB(
      Channel.value( file(params.sequences) ),
      Channel.value( file(params.seq2tax) ),
      Channel.value( file(params.taxonomy_list) )
    )
  }

  // ─── SUBCOMAND: COLLAPSE-TAXONOMY ────────────────────────────────────────
  else if( params.cmd == 'collapse-taxonomy' ) {
    if( !params.input_tsv )
      error "For collapse-taxonomy you must give --input_tsv"

    COLLAPSE(
      Channel.value( file(params.input_tsv) ),
      Channel.value( params.rank )
    )
  }

  // ─── SUBCOMAND: COMBINE-OUTPUTS ─────────────────────────────────────────
  else if( params.cmd == 'combine-outputs' ) {
    if( !params.input_dir )
      error "For combine-outputs you must give --input_dir"

    COMBINE(
      Channel.value( file(params.input_dir) ),
      Channel.value( params.rank )
    )
  }

  // ─── UNKNOWN COMMAND ───────────────────────────────────────────────────
  else {
    error "Unrecognized command: ${params.cmd}"
  }
}





