// main.nf – Nextflow version of the VI abundance pipeline
nextflow.enable.dsl = 2

// ─── INCLUDES (DSL2) ────────────────────────────────────────────────────────
include { FILTER_READS }  from './modules/abundance/filter_reads.nf'
include { ALIGN_READS }   from './modules/abundance/alignment.nf'
include { GET_CIGAR }     from './modules/abundance/cigar_probs.nf'
include { COMPUTE_LOGP }  from './modules/abundance/log_prob_rgs.nf'
include { RUN_VI }        from './modules/abundance/inference.nf'
include { WRITE_OUTPUT }  from './modules/abundance/write_output.nf'

include { BUILD_DB }      from './modules/build_database/main.nf'
include { COLLAPSE }      from './modules/collapse_taxonomy/main.nf'
include { COMBINE }       from './modules/combine_outputs/main.nf'

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

    // Step 2: Alignment
    fasta_file = file("${params.db}/species_taxid.fasta")
    sam_file = ALIGN_READS(
      named_reads,    // Tuple (sample_id, fastq.gz)
      fasta_file      // Input 
    )

    // Step 3: Get CIGAR probabilities
    cigar_info = GET_CIGAR(sam_file, params.threads)

    // ✅ Join SAM and CIGAR for sample_id
    sam_and_cigar = sam_file.join(cigar_info)

    // Step 4: Log-probabilidades calculate
    logp_data = COMPUTE_LOGP(sam_and_cigar)


    // Step 5: Variational inference algorithm
    freq_output = RUN_VI(logp_data)

    // Step 6: Get assigned/unassigned counts from JSON
    counts = logp_data.map { sample_id, file ->
      def json = file.text
      def parsed = new groovy.json.JsonSlurper().parseText(json)
      tuple(sample_id, parsed.assigned_count, parsed.unassigned_count)
    }

    assigned   = counts.map { it[1] }
    unassigned = counts.map { it[2] }


    // Step 7: Get final output
    WRITE_OUTPUT(
    freq_output,                    
    file(params.taxonomy_tsv),      
    assigned,                        
    unassigned                       
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





