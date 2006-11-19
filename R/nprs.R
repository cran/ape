### nprs.R (2003-07-11)
###
###    Nonparametric Rate Smoothing Method by Sanderson
###
### Copyright 2003 Gangolf Jobb and Korbinian Strimmer
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

setTree <-
  function(lowerNodes,upperNodes,edgeLengths,minEdgeLength,tipLabels)
   .C(
    "setTree",
    as.integer(lowerNodes),
    as.integer(upperNodes),
    as.double(edgeLengths),
    as.double(minEdgeLength),
    as.integer(length(edgeLengths)),
    as.character(tipLabels),
    as.integer(length(tipLabels)),
    result=integer(1),
    PACKAGE = "ape"
   )$result

getNFreeParams <-
  function()
   .C(
    "getNFreeParams",
    result=integer(1),
    PACKAGE = "ape"
   )$result

getNEdges <-
  function()
   .C(
    "getNEdges",
    result=integer(1),
    PACKAGE = "ape"
   )$result

getEdgeLengths <-
  function()
   .C(
    "getEdgeLengths",
    result=double(getNEdges()),
    PACKAGE = "ape"
   )$result

objFuncLogScale <-
  function(params,expo)
   .C(
    "objFuncLogScale",
    as.double(params),
    as.integer(expo),
    result=double(1),
    PACKAGE = "ape"
   )$result

getDurations <-
  function(params,scale)
   .C(
    "getDurations",
    as.double(params),
    as.double(scale),
    result=double(getNEdges()),
    PACKAGE = "ape"
   )$result

getRates <-
  function(params,scale)
   .C(
    "getRates",
    as.double(params),
    as.double(scale),
    result=double(getNEdges()),
    PACKAGE = "ape"
   )$result

getExternalParams <-
  function()
   .C(
    "getExternalParams",
    result=double(getNFreeParams()),
    PACKAGE = "ape"
   )$result

### private functions

prepareTree <- function(phy, minEdgeLength = 1e-06)
{
    len <- phy$edge.length
    if (length(len) > 2048) stop("Only 2048 branches in tree allowed!")
    low <- phy$edge[, 1] # edges in the tree
    upp <- phy$edge[, 2]
    setTree(low, upp, len, minEdgeLength, phy$tip.labels)
}

optimTree <- function(phy, expo = 2) # call prepareTree first
{
    dur <- rep(log(0.5), getNFreeParams() ) # start value
    objL <- function(d) objFuncLogScale(d, expo)
    opt <- optim(dur, objL, method = "BFGS")
    return(opt)
}

### this is just for testing purposes, to get the tree we are
### actually using when there are many small branch lengths
phylogram <- function(phy, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")

    ## added by EP for the new coding of "phylo" (2006-10-04):
    phy <- new2old.phylo(phy)
    ## End

    prepareTree(phy, ...)
    ##opt <- optimTree(phy, ...)

    newTree <- phy
    newTree$edge.length <- getEdgeLengths()

    ans <- newTree
    old2new.phylo(ans)
}

### public functions

chronogram <- function(phy, scale = 1, expo = 2, minEdgeLength = 1e-06)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")

    ## added by EP for the new coding of "phylo" (2006-10-04):
    phy <- new2old.phylo(phy)
    ## End

    prepareTree(phy, minEdgeLength = minEdgeLength)
    opt <- optimTree(phy, expo = expo)

    newTree <- phy
    newTree$edge.length <- getDurations(opt$par, scale)

    ans <- newTree
    old2new.phylo(ans)
}

ratogram <- function(phy, scale = 1, expo = 2, minEdgeLength = 1e-06)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")

    ## added by EP for the new coding of "phylo" (2006-10-04):
    phy <- new2old.phylo(phy)
    ## End

    prepareTree(phy, minEdgeLength = minEdgeLength)
    opt <- optimTree(phy, expo = expo)

    newTree <- phy
    newTree$edge.length <- getRates(opt$par, scale)

    ans <- newTree
    old2new.phylo(ans)
}

NPRS.criterion <- function(phy, chrono, expo = 2, minEdgeLength = 1e-06)
{
    if (!is.ultrametric(chrono))
      stop("tree \"chrono\" is not ultrametric (clock-like)")

    ## added by EP for the new coding of "phylo" (2006-10-04):
    phy <- new2old.phylo(phy)
    chrono <- new2old.phylo(chrono)
    ## End

    prepareTree(chrono, minEdgeLength = minEdgeLength)
    parms <- getExternalParams()
    prepareTree(phy, minEdgeLength = minEdgeLength)
    objFuncLogScale(parms, expo)
}
