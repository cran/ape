## pic.R (2009-05-10)

##   Phylogenetically Independent Contrasts

## Copyright 2002-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

pic <- function(x, phy, scaled = TRUE, var.contrasts = FALSE)
{
    if (!inherits(phy, "phylo"))
      stop("object 'phy' is not of class \"phylo\"")
    if (is.null(phy$edge.length))
      stop("your tree has no branch lengths: you may consider setting them equal to one, or using the function `compute.brlen'.")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
      stop("'phy' is not rooted and fully dichotomous")
    if (length(x) != nb.tip)
      stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x)))
      stop("the present method cannot (yet) be used directly with missing data: you may consider removing the species with missing data from your tree with the function `drop.tip'.")

    phy <- reorder(phy, "pruningwise")
    phenotype <- numeric(nb.tip + nb.node)

    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    } else {
        if (all(names(x) %in% phy$tip.label))
          phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning('the names of argument "x" and the tip labels of the tree did not match: the former were ignored in the analysis.')
        }
    }
    ## No need to copy the branch lengths: they are rescaled
    ## in the C code, so it's important to leave the default
    ## `DUP = TRUE' of .C.
    contr <- var.con <- numeric(nb.node)

    ans <- .C("pic", as.integer(nb.tip), as.integer(nb.node),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              as.double(phy$edge.length), as.double(phenotype),
              as.double(contr), as.double(var.con),
              as.integer(var.contrasts), as.integer(scaled),
              PACKAGE = "ape")

    ## The "old" R code:
    ##for (i in seq(from = 1, by = 2, length.out = nb.node)) {
    ##    j <- i + 1
    ##    anc <- phy$edge[i, 1]
    ##    des1 <- phy$edge[i, 2]
    ##    des2 <- phy$edge[j, 2]
    ##    sumbl <- bl[i] + bl[j]
    ##    ic <- anc - nb.tip
    ##    contr[ic] <- phenotype[des1] - phenotype[des2]
    ##    if (scaled) contr[ic] <- contr[ic]/sqrt(sumbl)
    ##    if (var.contrasts) var.con[ic] <- sumbl
    ##    phenotype[anc] <- (phenotype[des1]*bl[j] + phenotype[des2]*bl[i])/sumbl
    ##    k <- which(phy$edge[, 2] == anc)
    ##    bl[k] <- bl[k] + bl[i]*bl[j]/sumbl
    ##
    ##}
    contr <- ans[[7]]
    if (var.contrasts) {
        contr <- cbind(contr, ans[[8]])
        dimnames(contr) <- list(1:nb.node + nb.tip, c("contrasts", "variance"))
    } else names(contr) <- 1:nb.node + nb.tip
    contr
}
