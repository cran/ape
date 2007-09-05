## klastorin.R (2003-05-26)

##   Klastorin's (1982) classifification method, applied to
##   phylogenetic trees as suggested by Misawa and Tajima (2000)

## Copyright 2003 Gangolf Jobb

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

########### PRIVATE ##############

getMisawaTajima <- function()
  .C("getMisawaTajima", result = integer(nTips()),
     PACKAGE = "ape")$result

### functions to set and extract phylo tree ###

buildTreeFromPhylo <- function(tree) {
    lowerNodes <- tree$edge[,1]
    upperNodes <- tree$edge[,2]
    edgeLengths <- tree$edge.length
    tipLabels <- tree$tip.label
    .C("buildTreeFromPhylo", as.integer(lowerNodes),
       as.integer(upperNodes), as.double(edgeLengths),
       as.integer(length(edgeLengths)),
       as.character(tipLabels),
       as.integer(length(tipLabels)),
       result = integer(1), PACKAGE = "ape"
       )$result
}

destroyTree <- function()
  .C("destroyTree", result = integer(1),
     PACKAGE = "ape")$result

getError <- function()
  .C("getError", result = integer(1),
     PACKAGE = "ape")$result

nTips <- function()
  .C("nTips", result = integer(1),
     PACKAGE = "ape")$result

########### PUBLIC ##############

klastorin <- function(phy)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    ## added by EP for the new coding of "phylo" (2006-10-04):
    phy <- new2old.phylo(phy)
    ## End
    buildTreeFromPhylo(phy)
    if (getError() !=0) stop("Could not load \"phylo\" object")
    tmp <- getMisawaTajima()
    destroyTree()
    tmp
}
