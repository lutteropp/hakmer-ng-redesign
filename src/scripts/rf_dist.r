#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#run like this: Rscript --vanilla rf_dist.r tree1.txt tree2.txt
# or like this: Rscript --vanilla rf_dist.r tree1.txt tree2.txt 0.001

library(phangorn)
library(ape)

bipart <- function(x) 
{
  x <- reorder(x, "postorder")
  nTips <- as.integer(length(x$tip.label))
  .Call('_phangorn_bipartCPP', PACKAGE = 'phangorn', x$edge, nTips)
}


SHORTwise <- function (x, nTips, delete=FALSE) 
{
  v <- 1:nTips
  l <- lengths(x)
  lv <- floor(nTips/2)  
  for (i in seq_along(x)) { 
    if(l[i]>lv){
      y <- x[[i]]
      x[[i]] <- v[-y]
    }        
    if(l[i]==nTips/2){ 
      y <- x[[i]]
      if (y[1] != 1) 
        x[[i]] <- v[-y]
    }
  }
  if(any(l==nTips) && delete){
    x <- x[l!=nTips]
  }
  x
}

RF0 <- function(tree1, tree2=NULL, check.labels = TRUE, rooted=FALSE){   
  if(has.singles(tree1)) tree1 <- collapse.singles(tree1)
  if(has.singles(tree2)) tree2 <- collapse.singles(tree2)
  r1 <- is.rooted(tree1)
  r2 <- is.rooted(tree2)
  if(!rooted){
    if(r1) {
      tree1<-unroot(tree1)
      r1 <- FALSE
    }
    if(r2) {
      tree2<-unroot(tree2)
      r2 <- FALSE
    }
  }
  else{
    if(r1 != r2) {
      message("one tree is unrooted, unrooted both")
      tree1<-unroot(tree1)
      tree2<-unroot(tree2)
      r1 <- r2 <- FALSE
    }
  }  
  if (check.labels) {
    ind <- match(tree1$tip.label, tree2$tip.label)
    if (any(is.na(ind)) | length(tree1$tip.label) !=
        length(tree2$tip.label))
      stop("trees have different labels")
    tree2$tip.label <- tree2$tip.label[ind]
    #       tree2$edge[match(ind, tree2$edge[, 2]), 2] <- 1:length(ind)
    ind2 <- match(seq_along(ind), tree2$edge[, 2])
    tree2$edge[ind2, 2] <- order(ind)
  }
  #if(!is.binary(tree1) | !is.binary(tree2))message("Trees are not binary!")
  bp1 <- bipart(tree1)
  bp2 <- bipart(tree2)
  nTips <- length(tree1$tip.label)
  if(!rooted){
    bp1 <- SHORTwise(bp1, nTips)
    bp2 <- SHORTwise(bp2, nTips)    
  }
  RF <- sum(match(bp1, bp2, nomatch=0L)==0L) + sum(match(bp2, bp1, nomatch=0L)==0L)
  
  print("Nnode(tree1):")
  print(Nnode(tree1))
  
  print("bipart(tree1):")
  print(bp1)
  
  print("Nnode(tree2):")
  print(Nnode(tree2))
  
  print("bipart(tree2):")
  print(bp2)
  
  print("Plain RF-distance:")
  print(RF)
  RF <- RF / (Nnode(tree1) + Nnode(tree2) - 2)
  print("Normalized RF-distance:")
  print(RF)
  print("False positives:")
  print(sum(match(bp1, bp2, nomatch=0L)==0L))
  print("False negatives:")
  print(sum(match(bp2, bp1, nomatch=0L)==0L))
}

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

RF0(tree1, tree2)
