# Samples with no aligned read (F45)

Two FASTQ files for the `e2e` nf-test "Samples with no aligned read get an
all-Unassigned table beside a normal sample", which runs them beside the e2e
fixture's sample through a samplesheet. They are kept out of `tests/test_data/e2e/reads`,
which the `test` profile reads as a directory.

| File | Reads | What happens |
|---|---|---|
| `non16s.fastq.gz` | the 30 `non16S_*` reads of `tests/test_data/e2e/reads/e2e_mock.fastq.gz` (random sequence) | pass the length filter, align to nothing: 30 Unassigned |
| `too_short.fastq.gz` | the first 20 16S reads of the same file, cut to 400 bp | all fail `--min_length 500`; fastplong writes an empty file, minimap2 a header-only BAM: 0 reads |
