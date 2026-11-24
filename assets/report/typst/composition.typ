#import "lib.typ": *
#let s = sec("composition")
#let f = s.facts

= Composition
Read-level relative abundance of each sample's main species (#figref(s, "composition_species")) and genera (#figref(s, "composition_genus")): the union of each sample's top #f.top, at most #f.cap, the rest as Other#if f.heatmap [; #figref(s, "heatmap_species") shows the #f.heatmap_taxa species with the highest mean abundance in every sample]. #tabref(s, "prevalence") lists the species detected in the most samples.

#fig(s, "composition_species")
#fig(s, "composition_genus")
#fig(s, "heatmap_species")
#datatable(s, "prevalence")
