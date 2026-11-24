#import "lib.typ": *
#let s = sec("flagged")
#let f = s.facts

= Flagged species
A species is flagged in a sample when its shared_read_fraction is at least #f.cutoff (most of its expected reads are split with another reported species) or when it shares an identical 16S with other species of its genus (identical_16s_species). #f.n_flagged species are flagged in at least one sample: #f.n_shared by shared_read_fraction and #f.n_identical by an identical 16S (#figref(s, "flagged_species"), #tabref(s, "flagged")). Their abundance rests on reads that 16S does not assign to them alone.

#fig(s, "flagged_species")
#datatable(s, "flagged")
