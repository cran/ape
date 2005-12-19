### as.phylo.formula.R  (2005-12-10)
###
###    Conversion from Taxonomy Variables to Phylogenetic Trees
###
### Copyright 2005 Julien Dutheil <julien.dutheil@univ-montp2.fr>
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

as.phylo.formula <- function(x, data=parent.frame(), ...)
{
  # Testing formula syntax:
  err <- "Formula must be of the kind \"~A1/A2/.../An\"."
  if(length(x) != 2) stop(err)
  if(x[[1]] != "~") stop(err)
  f <- x[[2]]
  taxo <- list()
  while(length(f) == 3) {
    if(f[[1]] != "/") stop(err)
    if(!is.factor(data[[deparse(f[[3]])]])) stop(paste("Variable", deparse(f[[3]]), "must be a factor."))
    taxo[[deparse(f[[3]])]] <- data[[deparse(f[[3]])]]
    if(length(f) > 1) f <- f[[2]]
  }
  if(!is.factor(data[[deparse(f)]])) stop(paste("Variable", deparse(f), "must be a factor."))
  taxo[[deparse(f)]] <- data[[deparse(f)]]
  taxo.data <- as.data.frame(taxo)
  leaves.names <- as.character(taxo.data[,1])
  taxo.data[,1] <- 1:nrow(taxo.data)
  # Now builds the phylogeny:

  f.rec <- function(subtaxo) { # Recurrent utility function
    u <- ncol(subtaxo)
    levels <- unique(subtaxo[,u])
    if(u == 1) {
      if(length(levels) != nrow(subtaxo))
        warning("Error, leaves names are not unique.")
      return(as.character(subtaxo[,1]))
    }
    t <- character(length(levels))
    for(l in 1:length(levels)) {
      x <- f.rec(subtaxo[subtaxo[,u] == levels[l],][1:(u-1)])
      if(length(x) == 1) t[l] <- x
      else t[l] <- paste("(", paste(x, collapse=","), ")", sep="")
    }
    return(t)
  }
  string <- paste("(", paste(f.rec(taxo.data), collapse=","), ");", sep="")
  phy<-read.tree(text=string)
  phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
  return(phy)
}
