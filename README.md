# Scripts and workflows used to produce figures

## Prerequisites

- Nextflow
- Singularity

## General commands

Rhe file `run.sh` indicates the commands to run. 

## Inferring the tree

```
nextflow run gofasta.nf -c gofasta.config
```

It will produce the output files in the `results` folder, in particular:

- ba286.aligned.fasta_closest.txt: The name of the closest sequences to BA.2.86 queries
- ba286.aligned.fasta_closest_unique_sequences.fasta.gz: The sequences of the closest sequences to BA.2.86 queries
- ba286.aligned.fasta_closest_unique_sequences.fasta_masked.fasta.gz: The masked sequences of the closest sequences to BA.2.86 queries
- context.aligned.fasta_masked.fasta.gz: The masked sequences of the contextual samples
- gisaid_pangolin.tsv: Pangolin annotation of all the sequences
- bootaligns/*: Bootstrap alignments
- boottrees/*: Bootstrap trees
- align.treefile: Phylogenetic tree
- align_fbp.treefile: Phylogenetic tree with bootstrap supports

## Producing the figures

R
```
> source("figures/lineageStats.R")
> source("figures/plot_mutations.R")
```

