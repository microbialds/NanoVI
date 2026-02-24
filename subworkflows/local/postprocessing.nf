// subworkflows/local/postprocessing.nf — Collapse + combine subworkflows

include { COLLAPSE } from '../../modules/local/collapse_taxonomy'
include { COMBINE }  from '../../modules/local/combine_outputs'

workflow COLLAPSE_TAXONOMY {

    take:
    input_tsv   // path: abundance TSV
    rank        // val: taxonomic rank
    bin_dir     // path: bin/ directory

    main:
    COLLAPSE(input_tsv, rank, bin_dir)

    emit:
    collapsed = COLLAPSE.out.collapsed
}

workflow COMBINE_OUTPUTS {

    take:
    input_dir   // path: directory with abundance tables
    rank        // val: taxonomic rank
    bin_dir     // path: bin/ directory

    main:
    COMBINE(input_dir, rank, bin_dir)

    emit:
    combined = COMBINE.out.combined
}
