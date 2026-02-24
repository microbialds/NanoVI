// subworkflows/local/abundance.nf — Abundance subworkflow

include { FILTER_READS }  from '../../modules/local/filter_reads'
include { INDEX_DB }      from '../../modules/local/index_db'
include { ALIGN_READS }   from '../../modules/local/align_reads'
include { GET_CIGAR }     from '../../modules/local/cigar_probs'
include { COMPUTE_LOGP }  from '../../modules/local/log_prob_rgs'
include { RUN_VI }        from '../../modules/local/inference'
include { WRITE_OUTPUT }              from '../../modules/local/write_output'
include { GENERATE_PHYLOSEQ_TABLES }  from '../../modules/local/generate_phyloseq_tables'

workflow ABUNDANCE {

    take:
    input_reads     // channel: tuple(sample_id, fastq)
    fasta_file      // path: reference FASTA
    taxonomy_tsv    // path: taxonomy TSV

    main:
    bin_dir = file("${projectDir}/bin")

    // Step 1: Filter reads by length
    named_reads = FILTER_READS(input_reads).filtered

    // Step 2: Build minimap2 index (once)
    db_index = INDEX_DB(fasta_file)

    // Step 3: Align reads against reference
    sam_file = ALIGN_READS(named_reads, db_index)

    // Step 4: Compute CIGAR operation log-probabilities
    cigar_info = GET_CIGAR(sam_file, bin_dir)

    // Step 5: Join SAM and CIGAR by sample_id, compute log-probabilities
    sam_and_cigar = sam_file.join(cigar_info)
    logp_data = COMPUTE_LOGP(sam_and_cigar, bin_dir)

    // Step 6: Run variational inference
    freq_output = RUN_VI(logp_data, bin_dir)

    // Step 7: Join abundance with logp data and write final output
    vi_with_counts = freq_output.join(logp_data)
    WRITE_OUTPUT(vi_with_counts, taxonomy_tsv, bin_dir)

    // Step 8: Consolidate per-sample outputs into phyloseq-compatible tables
    collected_tsvs = WRITE_OUTPUT.out.result_tsv.collect()
    GENERATE_PHYLOSEQ_TABLES(collected_tsvs, bin_dir)

    emit:
    results        = WRITE_OUTPUT.out.result_tsv
    otu_abundance  = GENERATE_PHYLOSEQ_TABLES.out.otu_abundance
    otu_counts     = GENERATE_PHYLOSEQ_TABLES.out.otu_counts
    tax_table      = GENERATE_PHYLOSEQ_TABLES.out.tax_table
    sample_metadata = GENERATE_PHYLOSEQ_TABLES.out.sample_metadata
}
