#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(phangorn)
library(ape)

if (length(args) != 2 && length(args) != 3) {
  stop("Usage: Rscript --vanilla rf_dist.r TREE_1_PATH TREE_2_PATH [MINIMUM_BRANCH_LENGTH]\n", call.=FALSE)
}

print(args[1])
print(args[2])
minbr <- 0

tree1 <- read.tree(file=args[1])
tree2 <- read.tree(file=args[2])

if (length(args) == 3) {
  minbr <- args[3]
  tree1 <- di2multi(tree1, tol = minbr)
  tree2 <- di2multi(tree2, tol = minbr)
}

print(RF.dist(tree1, tree2, normalize = TRUE))

#run like this: Rscript --vanilla rf_dist.r tree1.txt tree2.txt
# or like this: Rscript --vanilla rf_dist.r tree1.txt tree2.txt 0.001
