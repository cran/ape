### klastorin.R  (2003-05-26)
###
###     Klastorin's (1982) classifification method, applied to 
###     phylogenetic trees as suggested by Misawa and Tajima (2000)
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
    buildTreeFromPhylo(phy) 
    if (getError() !=0) stop("Could not load \"phylo\" object")
    tmp <- getMisawaTajima()
    destroyTree()
    tmp
}
