// subworkflows/local/abundance.nf — Abundance subworkflow

include { FILTER_READS }  from '../../modules/local/filter_reads'
include { INDEX_DB }      from '../../modules/local/index_db'
include { ALIGN_READS }   from '../../modules/local/align_reads'
include { LIKELIHOOD }    from '../../modules/local/likelihood'
include { RUN_VI }        from '../../modules/local/inference'
include { IDENTICAL_16S }             from '../../modules/local/identical_16s'
include { WRITE_OUTPUT }              from '../../modules/local/write_output'
include { GENERATE_PHYLOSEQ_TABLES }  from '../../modules/local/generate_phyloseq_tables'

workflow ABUNDANCE {

    take:
    input_reads     // channel: tuple(sample_id, fastq)
    fasta_file      // path: reference FASTA (not used when db_index is given)
    taxonomy_tsv    // path: taxonomy TSV
    db_index        // path: prebuilt minimap2 index (--db_index), or null to build it
    identical_16s   // path: the database's identical_16s.tsv, or null to build it from the FASTA

    main:
    // Step 1: Filter reads by length
    named_reads = FILTER_READS(input_reads).filtered

    // A stub run's empty outputs get their own names, so a later real run never reuses them
    def stub = workflow.stubRun ? '.stub' : ''
    // The FASTA's MD5 names the cached index and identical-16S table (read once; skipped when neither is built)
    def fasta_md5 = (!db_index || !identical_16s) && fasta_file.exists() ? md5Prefix(fasta_file) : null

    // Step 2: minimap2 index: the prebuilt one, or built once per database, preset,
    // k-mer size and minimap2 version and kept in the index cache (F26)
    if (db_index) {
        mm2_index = channel.value(db_index)
        index_versions = channel.empty()
    } else {
        def index_name = "${fasta_file.baseName}.${fasta_md5}.${params.type}.k${params.kmer_size}${stub}"
        INDEX_DB(fasta_file, index_name)
        mm2_index = INDEX_DB.out.index
        index_versions = INDEX_DB.out.versions
    }

    // Species of one genus that share an identical 16S (D24): the database's table, or
    // built once from its FASTA and taxonomy and kept in the index cache. Without
    // either, the identical_16s_species column is left empty
    if (identical_16s) {
        identical_table = channel.value(identical_16s)
        identical_versions = channel.empty()
    } else if (fasta_md5) {
        // Named after the FASTA, the taxonomy and the script that groups them, so a
        // changed grouping rule builds a new table instead of reusing a cached one
        def script_md5 = md5Prefix(file("${projectDir}/bin/identical_16s.py"))
        def table_name = "${fasta_file.baseName}.${fasta_md5}.${md5Prefix(taxonomy_tsv)}.identical_16s-${script_md5}${stub}"
        IDENTICAL_16S(fasta_file, taxonomy_tsv, table_name)
        identical_table = IDENTICAL_16S.out.table
        identical_versions = IDENTICAL_16S.out.versions
    } else {
        log.warn "No identical_16s.tsv in --db and no species_taxid.fasta to build it from: " +
            "the identical_16s_species column is left empty"
        identical_table = channel.value([])
        identical_versions = channel.empty()
    }

    // Step 3: Align reads against reference
    bam_file = ALIGN_READS(named_reads, mm2_index).bam

    // Step 4: Error profile and per-read scores, from one pass over the alignments (F24)
    LIKELIHOOD(bam_file)

    // Step 5: Run variational inference on the per-read scores
    freq_output = RUN_VI(LIKELIHOOD.out.scores).abundance

    // Step 6: Join abundance with the read counts and write final output
    vi_with_counts = freq_output.join(LIKELIHOOD.out.counts)
    WRITE_OUTPUT(vi_with_counts, taxonomy_tsv, identical_table)

    // Step 7: Consolidate per-sample outputs into phyloseq-compatible tables
    collected_tsvs = WRITE_OUTPUT.out.result_tsv.collect()
    GENERATE_PHYLOSEQ_TABLES(collected_tsvs)

    // Software versions reported by every step (deduplicated in main.nf)
    versions = FILTER_READS.out.versions
        .mix(index_versions, identical_versions, ALIGN_READS.out.versions, LIKELIHOOD.out.versions,
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

// First 12 hex digits of the file's MD5, read in 1 MB blocks: identifies the
// database in the index file name, so two databases never share a cached index
def md5Prefix(path) {
    def digest = java.security.MessageDigest.getInstance('MD5')
    path.withInputStream { stream ->
        stream.eachByte(1 << 20) { buffer, n -> digest.update(buffer, 0, n) }
    }
    return digest.digest().encodeHex().toString().take(12)
}
