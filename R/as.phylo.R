### as.phylo.R  (2003-05-26)
###
###     Conversion between phylo and hclust trees
###
### Copyright 2003 Gangolf Jobb <gangolf@treefinder.de>
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

########### PRIVATE ##############

### functions to set and extract phylo tree ###

 buildTreeFromPhylo <-
  function(tree) {
   lowerNodes <- tree$edge[,1]
   upperNodes <- tree$edge[,2]
   edgeLengths <- tree$edge.length
   tipLabels <- tree$tip.label
   .C(
    "buildTreeFromPhylo",
    as.integer(lowerNodes),
    as.integer(upperNodes),
    as.double(edgeLengths),
    as.integer(length(edgeLengths)),
    as.character(tipLabels),
    as.integer(length(tipLabels)),
    result=integer(1),
    PACKAGE = "ape"
   )$result   
 }

 destroyTree <-
  function()
   .C(
    "destroyTree",
    result=integer(1),
    PACKAGE = "ape"
   )$result

 getError <-
  function()
   .C(
    "getError",
    result=integer(1),
    PACKAGE = "ape"
   )$result

 nTips <-
  function()
   .C(
    "nTips",
    result=integer(1),
    PACKAGE = "ape"
   )$result

 nNodes <-
  function()
   .C(
    "nNodes",
    result=integer(1),
    PACKAGE = "ape"
   )$result

 nEdges <-
  function()
   .C(
    "nEdges",
    result=integer(1),
    PACKAGE = "ape"
   )$result

 tipLabelsForPhylo <-
  function() {
   .C(
    "tipLabelsForPhylo",
    result=character(nTips()),
    PACKAGE = "ape"
   )$result
  }

 edgeLengthsForPhylo <-
  function() {
   .C(
    "edgeLengthsForPhylo",
    result=double(nEdges()),
    PACKAGE = "ape"
   )$result
  }

 lowerNodesForPhylo <-
  function() {
   .C(
    "lowerNodesForPhylo",
    result=integer(nEdges()),
    PACKAGE = "ape"
   )$result
  }

 upperNodesForPhylo <-
  function() {
   .C(
    "upperNodesForPhylo",
    result=integer(nEdges()),
    PACKAGE = "ape"
   )$result
  }

 getPhylo <-
  function() {
   edge <- matrix(nrow=nEdges(),ncol=2)
   edge[,1] <- lowerNodesForPhylo()
   edge[,2] <- upperNodesForPhylo()
   edge.length <- edgeLengthsForPhylo()
   tip.label <- tipLabelsForPhylo()
   mode(edge) <- "character"
   tree <- list(edge = edge,edge.length = edge.length,tip.label = tip.label)
   class(tree) <- "phylo"     
   return(tree)
  }


### functions to extract hclust tree ###

 buildTreeFromHclust <-
  function(tree) {
   leftNodes <- tree$merge[,1]
   rightNodes <- tree$merge[,2]
   nodeHeights <- tree$height
   tipLabels <- tree$labels
   tipOrder <- tree$order
   .C(
    "buildTreeFromHclust",
    as.integer(leftNodes),
    as.integer(rightNodes),
    as.double(nodeHeights),
    as.integer(length(nodeHeights)),
    as.character(tipLabels),
    as.integer(length(tipLabels)),
    result=integer(1),
    PACKAGE = "ape"
   )$result
 }

 leftNodesForHclust <-
  function() {
   .C(
    "leftNodesForHclust",
    result=integer(1+nEdges()-nTips()),
    PACKAGE = "ape"
   )$result
  }

 rightNodesForHclust <-
  function() {
   .C(
    "rightNodesForHclust",
    result=integer(1+nEdges()-nTips()),
    PACKAGE = "ape"
   )$result
  }

 nodeHeightsForHclust <-
  function() {
   .C(
    "nodeHeightsForHclust",
    result=double(1+nEdges()-nTips()),
    PACKAGE = "ape"
   )$result
  }

 getHclust <-
  function(mycall) {
   merge <- matrix(nrow=(1+nEdges()-nTips()),ncol=2)
   merge[,1] <- leftNodesForHclust()
   merge[,2] <- rightNodesForHclust()
   height <- nodeHeightsForHclust()
   labels <- tipLabelsForPhylo()
   order <- 1:length(labels)
   
   # sort according to heights
   #height.order <- order(height)
   #height <- height[height.order]
   #merge[,1] <- merge[height.order,1]
   #merge[,2] <- merge[height.order,2]
   
   #remap <- function(x, new.order)
   #{
     # remap all positive entries
   #  tmp <- x > 0 
    # x[tmp] <- new.order[x[tmp]]
    # return (x)
   #}
   #merge[,1] <- remap(merge[,1], height.order)
   #merge[,2] <- remap(merge[,2], height.order)
   
     
   storage.mode(merge) <- 'integer'
   storage.mode(order) <- 'integer'
      
   tree <- list(merge=merge,height=height, order=order,
                labels=labels,method='unknown',
		call=mycall
		)
   class(tree) <- "hclust"
   
   return(tree)
  }


########### PUBLIC ##############

# use method overloading 
# already defined in mva
#as.hclust <- function(phy) UseMethod("as.hclust")

# convert phylo object into hclust object
as.hclust.phylo <- function(x, ...)
{
  phy <- x
  rm(x)
  if (class(phy) != "phylo")
    stop("object is not of class \"phylo\"")

  if (is.ultrametric(phy) == FALSE) stop("object of class \"phylo\" not ultrametric") 
  if (is.binary.tree(phy) == FALSE) stop("object of class \"phylo\" not binary")
  
  buildTreeFromPhylo(phy)
  if (getError() !=0) stop("Could not load \"phylo\" object")
  
  hclass.out <- getHclust(match.call())
  
  destroyTree()
  
  return(hclass.out)
}

# use method overloading 
as.phylo <- function(hc) UseMethod("as.phylo")

# convert hclust object into phylo object
as.phylo.hclust <- function(hc)
{
  if (class(hc) != "hclust")
    stop("object \"hc\" is not of class \"hclust\"")

  buildTreeFromHclust(hc)  
  if (getError() !=0) stop("Could not load \"hclust\" object")
  
  phylo.out <- getPhylo()
  
  destroyTree()
  
  return(phylo.out)
}



