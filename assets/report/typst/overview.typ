#import "lib.typ": *
#let s = sec("overview")
#let f = s.facts

= Overview
This report covers #f.n_samples sample#if f.n_samples != 1 [s]#if f.raw_reads != "" [ with #f.raw_reads raw reads, of which #f.passed_reads (#f.pct_passed %) passed fastplong's length filter]#if f.has_counts [; #f.species_reads of the passed reads (#f.pct_species %) are on reported species, #f.below_reads below detection and #f.unassigned_reads unassigned]. #f.n_species species are reported in at least one sample (median #f.median_species per sample).

#if s.lists.attention.len() > 0 [
  == Attention
  #bullets(s.lists.attention)
]

== Read fate
#figref(s, "read_fate") and #tabref(s, "fate") follow each sample's reads from fastplong to the species table.

#fig(s, "read_fate")
#datatable(s, "fate")
#datatable(s, "metadata")
