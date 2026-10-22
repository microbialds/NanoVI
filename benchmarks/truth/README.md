# Truth tables of the real mock communities

Each table lists the species expected in a mock, under the name in its source (NCBI)
and under the GTDB r226 name that NanoVI's database uses (suffixes stripped; from
`benchmarks/names.py`, with the synonyms noted). `abundance` is the expected share of
16S reads, renormalised over the species that have a 16S gene; it is empty for
species known only to be present.

| File | Samples | Source |
|---|---|---|
| zymo_d6305_16s.tsv | ONT Zymo 2025.09 MSPlus (the 8 bacteria), SRR10391201, SRR23636353 | ZymoBIOMICS Microbial Community (DNA) Standard D6300/D6305 instruction manual, Table 1, "16S Only" column |
| zymo_msplus_spikes.tsv | ONT Zymo 2025.09 MSPlus | ONT release page (https://epi2me.nanoporetech.com/zymo_16s_2025.09/): four ATCC gDNA spike-ins, amounts not given |
| zymo_d6331_16s.tsv | ONT Zymo 2025.09 GutMock | ZymoBIOMICS Gut Microbiome Standard D6331 instruction manual, Table 1, "16S Only" column (five E. coli strains summed; the two yeasts have no 16S) |
| emu_gutmock_ont.tsv | SRR14326315 | Emu (Curry et al. 2022, Nature Methods) Supplementary Table 4, imputed ONT truth; Romboutsia hominis was absent from Emu's database |

The Emu gut mock's truth is itself an estimate (the authors imputed it from the reads
against a database of the strains in the mock), so it is a weaker reference than
Zymo's theoretical compositions.
