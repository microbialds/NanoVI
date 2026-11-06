#import "lib.typ": *
#let s = sec("reads")

= Per-read assignments
With \-\-read_assignments, each sample's file lists every read that passed fastplong with its table rows at probability 0.01 or more (D37). #tabref(s, "per_sample") counts the reads by their highest probability; #tabref(s, "pairs") lists the species pairs that share the most reads.

#datatable(s, "per_sample")
#datatable(s, "pairs")
