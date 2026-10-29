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
    // The FASTA's MD5 names the cached index and identical-16S table (skipped when neither is
    // built). It is read once and kept in the index cache, keyed by path, size and mtime (F50)
    def cache_dir = file(params.index_cache_dir ?: "${params.output_dir}/index_cache")
    def fasta_md5 = (!db_index || !identical_16s) && fasta_file.exists() ? cachedMd5Prefix(fasta_file, cache_dir) : null

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
        def table_name = "${fasta_file.baseName}.${fasta_md5}.${cachedMd5Prefix(taxonomy_tsv, cache_dir)}.identical_16s-${script_md5}${stub}"
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

    // A sample with no read aligned to the database (a blank, or every read filtered
    // out) gets a table with only the Below detection and Unassigned rows (F45)
    if (!workflow.stubRun) {
        LIKELIHOOD.out.counts.subscribe { sample_id, counts ->
            def n = new groovy.json.JsonSlurper().parseText(counts.text)
            if (!n.assigned_count)
                log.warn "Sample ${sample_id}: no read aligned to the database; the " +
                    "${n.unassigned_count} reads that passed the length filter are reported as Unassigned"
        }
    }

    // Step 5: Run variational inference on the per-read scores
    freq_output = RUN_VI(LIKELIHOOD.out.scores, taxonomy_tsv).abundance

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
    return md5Hex(path).take(12)
}

def md5Hex(path) {
    def digest = java.security.MessageDigest.getInstance('MD5')
    path.withInputStream { stream ->
        stream.eachByte(1 << 20) { buffer, n -> digest.update(buffer, 0, n) }
    }
    return digest.digest().encodeHex().toString()
}

// md5Prefix, reading the file only when <cache_dir>/md5_cache.tsv has no line for its
// path, size, modification time and, for a local file, inode and status-change time (F50):
// a GTDB-sized FASTA is read once, not at every launch, and a remote one is downloaded once.
// The names built from it do not change. The table is read again just before it is
// rewritten, so lines another launch added while this one read the file are kept, and it is
// rewritten through a temporary file, so a concurrent launch never reads half of it; when it
// cannot be read or written, the MD5 is still returned
def cachedMd5Prefix(path, cache_dir) {
    def table = cache_dir.resolve('md5_cache.tsv')
    def key = [path.toUriString(), path.size(), path.lastModified(), localFileId(path)].join('\t')
    def cached = readMd5Cache(table, "reading ${path} for its MD5")[key]
    if (cached ==~ /[0-9a-f]{32}/)
        return cached.take(12)
    def md5 = md5Hex(path)
    try {
        // A line for this path with another size or mtime can never match again
        def entries = readMd5Cache(table, "rewriting it").findAll { k, v -> !k.startsWith(path.toUriString() + '\t') }
        entries[key] = md5
        cache_dir.mkdirs()
        def tmp = cache_dir.resolve("md5_cache.tsv.${UUID.randomUUID()}.tmp")
        tmp.text = "path\tsize\tmtime_ms\tfile_id\tmd5\n" + entries.collect { k, v -> "${k}\t${v}\n" }.join('')
        java.nio.file.Files.move(tmp, table, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }
    catch (Exception e) {
        log.warn "Could not record the MD5 of ${path} in ${table}: ${e.message}"
    }
    return md5.take(12)
}

// The lines of md5_cache.tsv as key (path, size, mtime, file id) -> MD5; empty when the
// table is missing or cannot be read
def readMd5Cache(table, fallback) {
    def entries = [:]
    try {
        if (table.exists()) {
            table.readLines().drop(1).each { line ->
                def i = line.lastIndexOf('\t')
                if (i > 0) entries[line.substring(0, i)] = line.substring(i + 1)
            }
        }
    }
    catch (Exception e) {
        log.warn "Could not read ${table}: ${e.message}; ${fallback}"
    }
    return entries
}

// Inode and status-change time (ns) of a local file, '-' for a remote one. Users cannot set
// the status-change time, so a FASTA replaced with the same size and mtime (rsync -a,
// cp -p, tar x) still gets a new key
def localFileId(path) {
    try {
        def ino = java.nio.file.Files.getAttribute(path, 'unix:ino')
        def ctime = java.nio.file.Files.getAttribute(path, 'unix:ctime')
        return "${ino}:${ctime.to(java.util.concurrent.TimeUnit.NANOSECONDS)}"
    }
    catch (Exception e) {
        return '-'
    }
}
