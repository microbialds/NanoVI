// subworkflows/local/abundance.nf — Abundance subworkflow

include { FILTER_READS }  from '../../modules/local/filter_reads'
include { INDEX_DB }      from '../../modules/local/index_db'
include { ALIGN_READS }   from '../../modules/local/align_reads'
include { LIKELIHOOD }    from '../../modules/local/likelihood'
include { RUN_VI }        from '../../modules/local/inference'
include { WRITE_OUTPUT }              from '../../modules/local/write_output'
include { GENERATE_PHYLOSEQ_TABLES }  from '../../modules/local/generate_phyloseq_tables'

workflow ABUNDANCE {

    take:
    input_reads     // channel: tuple(sample_id, fastq)
    fasta_file      // path: reference FASTA
    taxonomy_tsv    // path: taxonomy TSV

    main:
    // Step 1: Filter reads by length
    named_reads = FILTER_READS(input_reads).filtered

    // Step 2: Build minimap2 index (once)
    db_index = INDEX_DB(fasta_file).index

    // Step 3: Align reads against reference
    sam_file = ALIGN_READS(named_reads, db_index).sam

    // Step 4: Error profile and per-read scores, from one pass over the alignments (F24)
    LIKELIHOOD(sam_file)

    // Step 5: Run variational inference on the per-read scores
    freq_output = RUN_VI(LIKELIHOOD.out.scores).abundance

    // Step 6: Join abundance with the read counts and write final output
    vi_with_counts = freq_output.join(LIKELIHOOD.out.counts)
    WRITE_OUTPUT(vi_with_counts, taxonomy_tsv)

    // Step 7: Consolidate per-sample outputs into phyloseq-compatible tables
    collected_tsvs = WRITE_OUTPUT.out.result_tsv.collect()
    GENERATE_PHYLOSEQ_TABLES(collected_tsvs)

    // Software versions reported by every step (deduplicated in main.nf)
    versions = FILTER_READS.out.versions
        .mix(INDEX_DB.out.versions, ALIGN_READS.out.versions, LIKELIHOOD.out.versions,
             RUN_VI.out.versions, WRITE_OUTPUT.out.versions,
             GENERATE_PHYLOSEQ_TABLES.out.versions)

    emit:
    results        = WRITE_OUTPUT.out.result_tsv
    otu_abundance  = GENERATE_PHYLOSEQ_TABLES.out.otu_abundance
    otu_counts     = GENERATE_PHYLOSEQ_TABLES.out.otu_counts
    tax_table      = GENERATE_PHYLOSEQ_TABLES.out.tax_table
    sample_metadata = GENERATE_PHYLOSEQ_TABLES.out.sample_metadata
    versions       = versions
}
