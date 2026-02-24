class SamplesheetParser {

    static List parseSamplesheet(samplesheet_path) {
        def rows = []
        def file = new File(samplesheet_path)

        if (!file.exists()) {
            throw new RuntimeException("Samplesheet not found: ${samplesheet_path}")
        }

        def lines = file.readLines()
        if (lines.size() < 2) {
            throw new RuntimeException("Samplesheet is empty or has no data rows: ${samplesheet_path}")
        }

        def header = lines[0].split(',').collect { it.trim().toLowerCase() }
        if (!header.contains('sample') || !header.contains('fastq')) {
            throw new RuntimeException("Samplesheet must have 'sample' and 'fastq' columns. Found: ${header}")
        }

        def sample_idx = header.indexOf('sample')
        def fastq_idx = header.indexOf('fastq')

        def seen_samples = [] as Set

        for (int i = 1; i < lines.size(); i++) {
            def line = lines[i].trim()
            if (line.isEmpty() || line.startsWith('#')) continue

            def fields = line.split(',').collect { it.trim() }
            def sample = fields[sample_idx]
            def fastq = fields[fastq_idx]

            // Validate sample name
            if (!(sample ==~ /^[a-zA-Z0-9_\-]+$/)) {
                throw new RuntimeException("Invalid sample name '${sample}' on line ${i+1}. Use alphanumeric, underscores, hyphens only.")
            }

            // Enforce unique sample IDs
            if (seen_samples.contains(sample)) {
                throw new RuntimeException("Duplicate sample ID '${sample}' found on line ${i+1}. Sample IDs must be unique.")
            }
            seen_samples.add(sample)

            // Require absolute paths
            if (!fastq.startsWith('/')) {
                throw new RuntimeException("FASTQ path must be absolute (start with '/'): ${fastq} (line ${i+1})")
            }

            // Validate FASTQ file exists
            def fq_file = new File(fastq)
            if (!fq_file.exists()) {
                throw new RuntimeException("FASTQ file not found: ${fastq} (line ${i+1})")
            }

            // Validate FASTQ extension
            if (!(fastq ==~ /.*\.(fastq|fq)(\.gz)?$/)) {
                throw new RuntimeException("Invalid FASTQ file extension: ${fastq} (line ${i+1}). Must end with .fastq, .fq, .fastq.gz, or .fq.gz")
            }

            rows.add([sample, fastq])
        }

        if (rows.isEmpty()) {
            throw new RuntimeException("No valid data rows in samplesheet: ${samplesheet_path}")
        }

        return rows
    }
}
