### aperates.R  (2003-05-05)
###
###     Nonparametric Rate Smoothing Method by Sanderson
###
### Copyright 2003 Gangolf Jobb <gangolf@treefinder.de> and
###                Korbinian Strimmer <strimmer@stat.uni-muenchen.de>
###
### This file is part of the `ape' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

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

    prepareTree(phy, ...)
    ##opt <- optimTree(phy, ...)
    
    newTree <- phy
    newTree$edge.length <- getEdgeLengths()
    
    return(newTree)
}

### public functions 

chronogram <- function(phy, scale = 1, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")

    prepareTree(phy, ...)
    opt <- optimTree(phy, ...)
    
    newTree <- phy
    newTree$edge.length <- getDurations(opt$par, scale)
        
    return(newTree)
}

ratogram <- function(phy, scale = 1, ...)
{
    if (class(phy) != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")

    prepareTree(phy, ...)
    opt <- optimTree(phy, ...)
    
    newTree <- phy
    newTree$edge.length <- getRates(opt$par, scale)
    
    return(newTree)
}

NPRS.criterion <- function(phy, chrono, expo = 2, ...)
{
    if (is.ultrametric(chrono) == FALSE)
      stop("tree \"chrono\" is not ultrametric (clock-like)")
    
    prepareTree(chrono, ...)
    parms <- getExternalParams() 
    prepareTree(phy, ...)
    objFuncLogScale(parms, expo) 
}
