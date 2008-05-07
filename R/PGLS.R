## PGLS.R (2008-05-08)

##   Phylogenetic Generalized Least Squares

## Copyright 2004 Julien Dutheil, and 2006-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

corBrownian <- function(value = 1, phy, form=~1)
{
  if (!("phylo" %in% class(phy))) stop("ERROR!!! Object \"phy\" is not of class \"phylo\"")
  attr(value, "formula") <- form
  attr(value, "fixed")   <- TRUE
  attr(value, "tree")    <- phy
  class(value) <- c("corBrownian", "corPhyl", "corStruct")
  return(value)
}

corMartins <- function(value, phy, form=~1, fixed=FALSE)
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

corGrafen <- function(value, phy, form=~1, fixed=FALSE)
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
  ## Added by EP 29 May 2006:
  if (nrow(data) != length(phy$tip.label))
    stop("number of observations and number of tips in the tree are not equal.")
  ## END
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
  mat <- vcv.phylo(tree, cor = corr)
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
  dist <- cophenetic.phylo(tree)
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
  mat <- vcv.phylo(tree, cor = corr)
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

### removed node.sons() and node.leafnumber()  (2006-10-12)

### changed by EP (2006-10-12):

compute.brlen <- function(phy, method = "Grafen", power = 1, ...)
{
    if (!"phylo" %in% class(phy))
      stop('object "phy" is not of class "phylo"')
    Ntip <- length(phy$tip.label)
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (is.numeric(method)) {
        phy$edge.length <- rep(method, length.out = Nedge)
        return(phy)
    }
    if (is.function(method)) {
        phy$edge.length <- method(Nedge, ...)
        return(phy)
    }
    if (is.character(method)) { # == "Grafen"
        tr <- reorder(phy, "pruningwise")
        xx <- .C("node_depth", as.integer(Ntip), as.integer(Nnode),
                 as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                 as.integer(Nedge), double(Ntip + Nnode),
                 DUP = FALSE, PACKAGE = "ape")[[6]] - 1
        m <- Ntip - 1
        phy$edge.length <-
          (xx[phy$edge[, 1]]/m)^power - (xx[phy$edge[, 2]]/m)^power
        return(phy)
    }
}

## by EP:

corPagel <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (value < 0 || value > 1)
        stop("the value of lambda must be between 0 and 1.")
    if (!("phylo" %in% class(phy)))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corPagel", "corPhyl", "corStruct")
    value
}

corMatrix.corPagel <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!any(attr(object, "index")))
        stop("object has not been initialized")
    mat <- vcv.phylo(attr(object, "tree"), cor = corr)
    index <- attr(object, "index")
    mat <- mat[index, index]
    tmp <- diag(mat)
    mat <- object[1]*mat
    diag(mat) <- tmp
    mat
}

coef.corPagel <- function(object, unconstrained = TRUE, ...)
{
    if (unconstrained) {
        if (attr(object, "fixed")) return(numeric(0))
        else return(object[1])
    }
    aux <- object[1]
    names(aux) <- "lambda"
    aux
}

corBlomberg <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (value <= 0)
        stop("the value of g must be greater than 0.")
    if (!("phylo" %in% class(phy)))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corBlomberg", "corPhyl", "corStruct")
    value
}

corMatrix.corBlomberg <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    index <- attr(object, "index")
    if (is.null(index))
        stop("object has not been initialized")
    if (object[1] <= 0)
        stop("the optimization has reached a value <= 0 for parameter 'g':
probably need to set 'fixed = TRUE' in corBlomberg().")
    phy <- attr(object, "tree")
    d <- (dist.nodes(phy)[length(phy$tip.label) + 1, ])^(1/object[1])
    phy$edge.length <- d[phy$edge[, 2]] - d[phy$edge[, 1]]
    mat <- vcv.phylo(phy, cor = corr)
    mat[index, index]
}

coef.corBlomberg <- function(object, unconstrained = TRUE, ...)
{
    if (unconstrained) {
        if (attr(object, "fixed")) return(numeric(0))
        else return(object[1])
    }
    aux <- object[1]
    names(aux) <- "g"
    aux
}
