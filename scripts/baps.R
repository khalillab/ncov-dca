#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("alignment",
                    help="MSA of isolates")
parser$add_argument("output",
                    help="BAPS clusters output")

args <- parser$parse_args()

suppressPackageStartupMessages(library("fastbaps"))
suppressPackageStartupMessages(library("ape"))

sparse.data <- import_fasta_sparse_nt(args$alignment)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
baps.hc <- fast_baps(sparse.data)
clusters <- best_baps_partition(sparse.data, as.phylo(baps.hc))

# save clusters
write.csv(as.data.frame(clusters),
	  file=args$output,
	  quote=FALSE)
