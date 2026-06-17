#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: plot_tree.R input_tree.nwk output.png")
}

treefile <- args[1]
outfile  <- args[2]

suppressPackageStartupMessages(library(ape))

tree <- read.tree(treefile)

png(outfile, width = 2000, height = 2000, res = 300)

# tree style options:
# “phylogram”, “cladogram”, “fan”, “unrooted”, “radial”, “tidy”

png(outfile, width = 2000, height = 2000, res = 300)
plot.phylo(
  tree,
  type = "unrooted",
  show.tip.label = FALSE,
  no.margin = TRUE
  # cex = 0.4,
)
dev.off()
