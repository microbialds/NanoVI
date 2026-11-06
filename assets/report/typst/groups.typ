#import "lib.typ": *
#let s = sec("groups")
#let f = s.facts

= Ambiguity groups
With \-\-ambiguity_groups, reported species of one family that share their reads form a group (D39). The #f.n_samples sample#if f.n_samples != 1 [s have] else [ has] #f.n_groups groups of two or more members, #f.n_sets distinct member sets (#figref(s, "groups_heatmap"), #tabref(s, "groups")). A group's abundance is more robust than the split between its members.

#fig(s, "groups_heatmap")
#datatable(s, "groups")
