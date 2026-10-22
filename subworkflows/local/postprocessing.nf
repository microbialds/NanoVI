// subworkflows/local/postprocessing.nf — Collapse + combine subworkflows

include { COLLAPSE } from '../../modules/local/collapse_taxonomy'
include { COMBINE }  from '../../modules/local/combine_outputs'

workflow COLLAPSE_TAXONOMY {

    take:
    input_tsv   // path: abundance TSV
    rank        // val: taxonomic rank

    main:
    COLLAPSE(input_tsv, rank)

    emit:
    collapsed = COLLAPSE.out.collapsed
    versions  = COLLAPSE.out.versions
}

workflow COMBINE_OUTPUTS {

    take:
    input_dir   // path: directory with abundance tables
    rank        // val: taxonomic rank

    main:
    COMBINE(input_dir, rank)

    emit:
    combined = COMBINE.out.combined
    versions = COMBINE.out.versions
}
