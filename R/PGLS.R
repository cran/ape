### PGLS.R  (2004-10-22)
###
###     Phylogenetic Generalized Least Squares
###
### Copyright 2004 Julien Dutheil <julien.dutheil@univ-montp2.fr>
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

# For debugging:
#library(ape)
#library(nlme)
#cat("((((Homo:0.21,Pongo:0.21):0.28,","Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",file = "ex.tre", sep = "\n")
#tree.primates <- read.tree("ex.tre")
#X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
#Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
#unlink("ex.tre") # delete the file "ex.tre"

corBrownian <- function(value = 1, phy, form=~1)
{
  if (!("phylo" %in% class(phy))) stop("ERROR!!! Object \"phy\" is not of class \"phylo\"")
  attr(value, "formula") <- form
  attr(value, "fixed")   <- TRUE
  attr(value, "tree")    <- phy
  class(value) <- c("corBrownian", "corPhyl", "corStruct")
  return(value)
}

corMartins <- function(value = numeric(0), phy, form=~1, fixed=FALSE)
{
  if(length(value) > 1) stop("ERROR!!! Only one parameter is allowed in corPGLS structure.")
  if(value < 0) stop("ERROR!!! Parameter alpha must be positive.")
  if (!("phylo" %in% class(phy))) stop("ERROR!!! Object \"phy\" is not of class \"phylo\"")
  attr(value, "formula") <- form
  attr(value, "fixed")   <- fixed
  attr(value, "tree")    <- phy
  class(value) <- c("corMartins", "corPhyl", "corStruct")
  return(value)
}

corGrafen <- function(value = numeric(0), phy, form=~1, fixed=FALSE)
{
  if(length(value) > 1) stop("ERROR!!! Only one parameter is allowed in corGrafen structure.")
  if(value < 0) stop("ERROR!!! Parameter rho must be positive.")
  value <- log(value) # Optimization under constraint, use exponential transform.
  if (!("phylo" %in% class(phy))) stop("ERROR!!! Object \"phy\" is not of class \"phylo\"")
  attr(value, "formula") <- form
  attr(value, "fixed")   <- fixed
  attr(value, "tree")    <- phy
  class(value) <- c("corGrafen", "corPhyl", "corStruct")
  return(value)
}

Initialize.corPhyl <- function(object, data, ...)
{
  # The same as in Initialize corStruct:
  form <- formula(object)
  ## Obtaining the group information, if any
  if(!is.null(getGroupsFormula(form))) {
    attr(object, "groups") <- getGroups(object, form, data = data)
    attr(object, "Dim")    <- Dim(object, attr(object, "groups"))
  } else { # no groups
    attr(object, "Dim")    <- Dim(object, as.factor(rep(1, nrow(data))))
  }
  ## Obtaining the covariate(s)
  attr(object, "covariate") <- getCovariate(object, data = data)

  # Specific to corPhyl:
  phy <- attr(object, "tree")
  if (is.null(data)) 
    data <- parent.frame()
  if(is.null(rownames(data))) {
    warning("No row names supplied in dataframe, data taken to be in the same order as in tree.")
    attr(object, "index") <- 1:dim(data)[1]
  } else {
    index <- match(rownames(data), phy$tip.label)
    if(any(is.na(index))) { 
      warning("Row names in dataframe do not match tree tip names. data taken to be in the same order as in tree.")
      attr(object, "index") <- 1:dim(data)[1]
    } else {
      attr(object, "index") <- index
    }
  }
  return(object)
}

corMatrix.corBrownian <- function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
  if (!("corBrownian" %in% class(object))) stop("ERROR!!! Object is not of class \"corBrownian\".")
  if(!any(attr(object, "index"))) stop("ERROR!!! object have not been initialized.")
  tree <- attr(object, "tree")
  mat <- vcv.phylo(tree, corr)
  n <- dim(mat)[1]
  # reorder matrix:
  matr <- matrix(nrow=n, ncol=n)
  index <- attr(object, "index")
  for(i in 1:n)
    for(j in i:n)
      matr[i,j] <- matr[j,i] <- mat[index[i], index[j]]
  return(matr)
}

corMatrix.corMartins <- function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
  if (!("corMartins" %in% class(object))) stop("ERROR!!! Object is not of class \"corMartins\".")
  if(!any(attr(object, "index"))) stop("ERROR!!! object have not been initialized.")
  tree <- attr(object, "tree")
  dist <- dist.phylo(tree)  
  mat <- exp(-object[1] * dist)
  if(corr) mat <- cov2cor(mat)
  n <- dim(mat)[1]
  # reorder matrix:
  matr <- matrix(nrow=n, ncol=n)
  index <- attr(object, "index")
  for(i in 1:n)
    for(j in i:n)
      matr[i,j] <- matr[j,i] <- mat[index[i], index[j]]
  return(matr)
}

corMatrix.corGrafen <- function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
  if (!("corGrafen" %in% class(object))) stop("ERROR!!! Object is not of class \"corGrafen\".")
  if(!any(attr(object, "index"))) stop("ERROR!!! object have not been initialized.")
  tree <- compute.brlen(attr(object, "tree"), method = "Grafen", power = exp(object[1]))
  mat <- vcv.phylo(tree, corr)
  n <- dim(mat)[1]
  # reorder matrix:
  matr <- matrix(nrow=n, ncol=n)
  index <- attr(object, "index")
  for(i in 1:n)
    for(j in i:n)
      matr[i,j] <- matr[j,i] <- mat[index[i], index[j]]
  return(matr)
}

coef.corBrownian <- function(object, unconstrained = TRUE, ...)
{
  if (!("corBrownian" %in% class(object))) stop("ERROR!!! Object is not of class \"corBrownian\".")
  return(numeric(0))  
}

coef.corMartins <- function(object, unconstrained = TRUE, ...)
{
  if (!("corMartins" %in% class(object))) stop("ERROR!!! Object is not of class \"corMartins\".")
  if(unconstrained) {
    if(attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux <- as.vector(object)
  names(aux) <- "alpha" 
  return(aux)
}

coef.corGrafen <- function(object, unconstrained = TRUE, ...)
{
  if (!("corGrafen" %in% class(object))) stop("ERROR!!! Object is not of class \"corGrafen\".")
  if(unconstrained) {
    if(attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }
  aux <- exp(as.vector(object))
  names(aux) <- "rho" 
  return(aux)
}

# Use Grafen's branch lengths:

node.sons <- function (phy, node)
{
    if (!("phylo" %in% class(phy)))
        stop("Object \"phy\" is not of class \"phylo\"")
    phy$edge[which(phy$edge[, 1] == node), 2]
}

node.leafnumber <- function(phy, node)
{
  if (!("phylo" %in% class(phy))) stop("Object \"phy\" is not of class \"phylo\"")
  if(as.numeric(node) > 0) return(1)
  else {
    number <- 0
    sons <- node.sons(phy, node)
    for(i in sons) number <- number + node.leafnumber(phy, i)
    return(number)
  }
}

compute.brlen <- function(phy, method="Grafen", power=1)
{
  if (!("phylo" %in% class(phy))) stop("Object \"phy\" is not of class \"phylo\"")
  E <- phy$edge
  n <- dim(E)[1]
  n.leaves <- length(phy$tip.label)
  for(i in 1:n) {
    bottom <- ((node.leafnumber(phy, E[i,1]) - 1)/(n.leaves - 1))^power
    top    <- ifelse(as.numeric(E[i, 2]) > 0, 0, ((node.leafnumber(phy, E[i, 2]) - 1)/(n.leaves - 1))^power)
    phy$edge.length[i] <- bottom - top 
  }
  #Now scale the tree:
  return(phy)
}

#m1 <- gls(Y~X, correlation=corBrownian(tree.primates))
#m2 <- gls(Y~X, correlation=corMartins(1, tree.primates))
#m3 <- gls(Y~X, correlation=corGrafen(1, tree.primates))

