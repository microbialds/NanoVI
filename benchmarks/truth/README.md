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

The Zymo standards' *Bacillus subtilis* is *Bacillus spizizenii* in GTDB r226: Zymo's
reference genome of the strain (ZymoBIOMICS.STD.refseq.v2) aligns at 99.97% identity
over its whole length to the *B. spizizenii* type strain ATCC 6633 (GCF_031316495.1,
GTDB's representative) and at 87.8% over 94% to *B. subtilis* 168 (GCF_000009045.1)
(minimap2 2.24 `-cx asm20`, alignments of MAPQ >= 20); five of its ten 16S copies are
identical only to *B. spizizenii* sequences of the database and the other five to
sequences shared with *B. subtilis* and relatives. The truth tables use *B. spizizenii*.

The Emu gut mock's truth is itself an estimate (the authors imputed it from the reads
against a database of the strains in the mock), so it is a weaker reference than
Zymo's theoretical compositions.
