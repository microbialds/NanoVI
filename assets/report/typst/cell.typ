#import "lib.typ": *
#let s = sec("cell")
#let f = s.facts

= Cell-level abundance
With --cell_abundance, each species' read-level abundance is divided by its 16S copy number and renormalised (D40)#if f.explanation != "" [: #f.explanation] else [.] #tabref(s, "left_out") lists what each sample's table leaves out; #figref(s, "cell_vs_read") compares the two abundances.

#datatable(s, "left_out")
#fig(s, "cell_vs_read")
#datatable(s, "cell")
