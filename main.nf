// main.nf - Nextflow version of the VI abundance pipeline

nextflow.enable.dsl = 2

include {FILTER_READS}  from "./modules/abundance/filter_reads.nf"
include {ALIGN_READS}   from "./modules/abundance/alignment.nf"
include {GET_CIGAR}     from "./modules/abundance/cigar_probs.nf"
include {COMPUTE_LOGP}  from "./modules/abundance/log_prob_rgs.nf"
include {RUN_VI}        from "./modules/abundance/inference.nf"
include {WRITE_OUTPUT}  from "./modules/abundance/write_output.nf"

workflow {

  if (params.cmd == 'abundance') {

    // Leer archivos FASTQ
    input_reads = Channel
      .fromPath("${params.input}/**", checkIfExists: true)
      .filter { it.name.endsWith('.fastq') || it.name.endsWith('.fq') || it.name.endsWith('.fastq.gz') || it.name.endsWith('.fq.gz') }
      .ifEmpty { error "❌ No FASTQ files found in ${params.input}" }

    // Step 1: Filtrado de lecturas (tupla: sample_id, fastq)
    named_reads = FILTER_READS(input_reads).filtered

    // Step 2: Alineamiento con la base como input explícito
    fasta_file = file("${params.db}/species_taxid.fasta")
    sam_file = ALIGN_READS(
      named_reads,    // Tupla (sample_id, fastq.gz)
      fasta_file      // Input como path explícito para evitar errores en Docker
    )

    // Step 3: Obtener probabilidades CIGAR
    cigar_info = GET_CIGAR(sam_file, params.threads)

    // Step 4: Calcular log-probabilidades
    logp_data = COMPUTE_LOGP(sam_file, cigar_info)

    // Step 5: Inferencia variacional
    freq_output = RUN_VI(logp_data)

    // Step 6: Obtener assigned/unassigned counts desde el JSON
    counts = logp_data.map { file ->
      def json = file.text
      def parsed = new groovy.json.JsonSlurper().parseText(json)
      tuple(parsed.assigned_count, parsed.unassigned_count)
    }

    assigned   = counts.map { it[0] }
    unassigned = counts.map { it[1] }

    // Step 7: Generar salida final
    WRITE_OUTPUT(
      abundance_json   = freq_output,
      taxonomy_tsv     = file(params.taxonomy_tsv),
      assigned_count   = assigned,
      unassigned_count = unassigned,
      output_prefix    = "${params.output_dir}/abundance_output"
    )
  }

  // Otros subcomandos
  else if (params.cmd == 'build-database') {
    include BUILD_DB from './modules/build_database/main.nf'
    BUILD_DB()
  }

  else if (params.cmd == 'collapse-taxonomy') {
    include COLLAPSE from './modules/collapse_taxonomy/main.nf'
    COLLAPSE()
  }

  else if (params.cmd == 'combine-outputs') {
    include COMBINE from './modules/combine_outputs/main.nf'
    COMBINE()
  }

  else {
    error "Comando no reconocido: ${params.cmd}"
  }
}







