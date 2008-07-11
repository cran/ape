## mlphylo.R (2008-06-18)

##   Estimating Phylogenies by Maximum Likelihood

## Copyright 2006-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

logLik.phylo <- function(object, ...) attr(object, "loglik")

deviance.phylo <- function(object, ...) -2*attr(object, "loglik")

AIC.phylo <- function(object, ..., k = 2)
{
    np <- length(object$edge.length) +
        length(attr(object, "rates")) +
            length(attr(object, "alpha")) +
                length(attr(object, "invar")) +
                    length(attr(object, "xi"))
    if (!attr(object, "model") %in% c("JC69", "F81"))
        np <- np + 3
    -2*attr(object, "loglik") + k*np
}

.subst.model <- structure(c(0, 1, 0, 1, 1, 1, 2, 5),
   names = c("JC69", "K80", "F81", "F84",
   "HKY85", "T92", "TN93", "GTR"))

mlphylo <-
    function(x, phy, model = DNAmodel(), search.tree = FALSE,
             quiet = FALSE, value = NULL, fixed = FALSE)
{
    ## not yet generic....
    if (class(x) != "DNAbin") stop("DNA sequences not in binary format")
    if (!is.binary.tree(phy))
        stop("the initial tree must be dichotomous.")
    if (!quiet && is.rooted(phy)) {
        warning("the initial tree is rooted: it will be unrooted.")
        phy <- unroot(phy)
    }
    if (is.null(phy$edge.length))
      stop("the initial tree must have branch lengths.")
    if (any(phy$edge.length > 1))
      stop("some branch lengths are greater than one.")
    phy <- reorder(phy, "pruningwise")
    if (!quiet) cat("Preparing the sequences...\n")
    if (is.list(x)) x <- as.matrix(x)
    if (is.null(rownames(x)))
        stop("DNA sequences have no names") # safe...
    if (!all(names(x) %in% phy$tip.label))
        stop("the names of the DNA sequences and the tip labels
of the tree do not match") # safe here also
    x <- x[phy$tip.label, ]
    Y <- prepareDNA(x, model)
    BF <- if (Y$model %in% 1:2) rep(0.25, 4) else base.freq(x)
    S <- length(Y$weight)
    npart <- dim(Y$partition)[2] # the number of overall partitions
    ## in case of negative branch lengths:
    phy$edge.length <- abs(phy$edge.length)
    nb.tip <- length(phy$tip.label)
    para <- if (Y$npara) rep(1, Y$npara) else 0
    alpha <- if (Y$nalpha) rep(.5, Y$nalpha) else 0
    invar <- if (Y$ninvar) rep(0.5, Y$ninvar) else 0

    if (!is.null(value)) {
        if (para && !is.null(value$rates))
            para <- value$rates[1:Y$npara]
        if (alpha && !is.null(value$alpha))
            alpha <- value$alpha[1:Y$nalpha]
        if (invar && !is.null(value$invar))
            invar <- value$invar[1:Y$ninvar]
    }

    loglik <- 0
    if (!quiet) cat("Fitting in progress... ")
    res <<- res <- .C("mlphylo_DNAmodel", as.integer(nb.tip), as.integer(S),
              as.raw(Y$SEQ), as.double(Y$ANC), as.double(Y$w),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              as.double(phy$edge.length), as.integer(npart),
              as.integer(Y$partition), as.integer(Y$model),
              as.double(Y$xi), as.double(para), as.integer(Y$npara),
              as.double(alpha), as.integer(Y$nalpha),
              as.integer(Y$ncat), as.double(invar), as.integer(Y$ninvar),
              as.double(BF), as.integer(search.tree), as.integer(fixed),
              as.double(loglik), NAOK = TRUE, PACKAGE = "mlphylo")
    if (!quiet) cat("DONE!\n")
    phy$edge.length = res[[8]]
    attr(phy, "loglik") <- res[[23]]
    attr(phy, "npart") <- npart
    attr(phy, "model") <- names(Y$npara)
    if (para) attr(phy, "rates") <- res[[13]]
    if (alpha) attr(phy, "alpha") <- res[[15]]
    if (invar) attr(phy, "invar") <- res[[18]]
    if (npart > 1) attr(phy, "xi") <- res[[12]]
    phy
}

DNAmodel <- function(model = "K80", partition = 1,
         ncat.isv = 1, invar = FALSE,
         equal.isv = TRUE, equal.invar = 1)
{
    if (ncat.isv > 10)
        stop("number of categories for inter-site variation cannot exceed 10")
    structure(list(model = model, partition = partition,
                   ncat.isv = ncat.isv, invar = invar,
                   equal.isv = equal.isv, equal.invar = equal.invar),
              class = "DNAmodel")
}

prepareDNA <- function(X, DNAmodel)
{
    L <- dim(X)[2] # already converted as a matrix in mlphylo()

    npart <- length(unique(DNAmodel$partition))

    ## find which substitution model:
    mo <- which(names(.subst.model) == DNAmodel$model)
    npara <- .subst.model[mo] # keeps the 'names'

    ## inter-sites variation:
    nalpha <- as.numeric(DNAmodel$ncat.isv > 1)
    if (!DNAmodel$equal.isv) nalpha <- npart * nalpha

    ## proportion of invariants:
    ninvar <- as.numeric(DNAmodel$invar)
    if (!DNAmodel$equal.invar) ninvar <- npart * ninvar

    SEQ <- weight <- part <- NULL

    ## For each partition...
    for (i in 1:npart) {
        ## extracts the sites in this partition:
        M <- X[, DNAmodel$partition == i, drop = FALSE]
        ## convert each column as a character string:
        M <- apply(M, 2, rawToChar)
        ## get their frequencies:
        w <- table(M)
        ## convert back to raw the unique(M):
        M <- sapply(dimnames(w)[[1]], charToRaw)
        ## remove useless attributes:
        colnames(M) <- dimnames(w) <- NULL
        w <- unclass(w)
        ## bind everything:
        SEQ <- cbind(SEQ, M)
        weight <- c(weight, w)
        part <- c(part, length(w)) # the length of each partition
    }

    class(SEQ) <- "DNAbin"
    ANC <- array(0, c(nrow(SEQ) - 2, ncol(SEQ), 4))

    ## 'partition' gives the start and end of each partition:
    partition <- matrix(1, 2, npart)
    partition[2, ] <- cumsum(part)
    if (npart > 1) {
        partition[1, 2:npart] <- partition[2, 1:(npart - 1)] + 1
        partition[2, npart] <- length(weight)
        xi <- rep(1, npart - 1)
    } else xi <- 0
    list(SEQ = SEQ, ANC = ANC, weight = weight, partition = partition,
         model = mo, xi = xi, npara = npara, nalpha = nalpha,
         ncat = DNAmodel$ncat.isv, ninvar = ninvar)
}
