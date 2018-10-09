#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(phangorn)

# test if there are at least two arguments: if not, return an error
if (length(args) != 2) {
  stop("Exactly two arguments must be supplied (tree input files).\n", call.=FALSE)
}

print(args[1])
print(args[2])

tree1 <- read.tree(file=args[1])
tree2 <- read.tree(file=args[2])

print(RF.dist(tree1, tree2, normalize = TRUE))

#run like this: Rscript --vanilla rf_dist.r tree1.txt tree2.txt
