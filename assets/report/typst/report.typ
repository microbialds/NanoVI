#import "lib.typ": *
#show: setup

#block(below: 4pt, text(19pt, weight: "bold", d.run.title))
#text(9.5pt, fill: muted)[Generated #d.run.generated#if d.run.nanovi_version != "" [ from a NanoVI #d.run.nanovi_version run] of #d.run.n_samples sample#if d.run.n_samples != 1 [s]#if d.run.nextflow_version != "" [ (Nextflow #d.run.nextflow_version)].]
#v(6pt)
#outline(depth: 2, indent: 1em)

// one file per section present, in the order Python gives
#for s in d.sections { include (s.id + ".typ") }
