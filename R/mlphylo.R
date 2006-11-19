### mlphylo.R (2006-07-05)
###
###     Estimating Phylogenies by Maximum Likelihood
###
### Copyright 2006 Emmanuel Paradis
###
### This file is part of the R-package `ape'.
### See the file ../COPYING for licensing issues.

logLik.phylo <- function(object, ...) attr(object, "loglik")

deviance.phylo <- function(object, ...) -2*attr(object, "loglik")

AIC.phylo <- function(object, ..., k = 2)
{
    np <- length(object$edge.length)
    PARA <- attr(object, "para")
    npart <- length(PARA)
    unbalancedBF <- FALSE
    for (i in 1:npart) {
        if (!PARA[[i]][[1]] %in% c("JC69", "F81"))
          np <- np + length(PARA[[i]][[2]])
        if (!PARA[[i]][[1]] %in% c("JC69", "K80"))
          unbalancedBF <- TRUE
    }
    if (unbalancedBF) np <- np + 3
    if (npart > 1) np <- np + npart - 1
    -2*attr(object, "loglik") + k * np
}

.subst.model <- c("JC69", "K80", "F81", "F84",
                  "HKY85", "T92", "TN93", "GTR")

mlphylo <- function(model = DNAmodel(), x, phy, search.tree = FALSE,
                    quiet = FALSE)
{
    if (!is.binary.tree(phy)) stop("the initial tree must be dichotomous.")
    if (is.rooted(phy))
      warning("the initial tree is rooted: it will be unrooted.")
    if (is.null(phy$edge.length))
      stop("the initial tree must have branch lengths.")
    if (!quiet) cat("Preparing the sequences...\n")
    BF <- base.freq(x)
    ## <FIXME> Will need to do the usual checks of names...
    x <- if (is.list(x)) x[phy$tip.label] else x[phy$tip.label, ]
    ## </FIXME>
    Y <- prepare.dna(x, model)
    S <- length(Y$weight)
    npart <- dim(Y$partition)[2] # the number of overall partitions

    ## prepare the tree:
    if (!quiet) cat("Preparing the tree...\n")
    if (!is.rooted(phy)) phy <- multi2di(phy, random = FALSE)

    ## in case NJ returns negative branch lengths:
    if (any(phy$edge.length < 0)) phy$edge.length <- abs(phy$edge.length)

    match <- as.matching(phy)
    nb.tip <- length(match$edge.length) / 2 + 1
    ## unroot the tree:
    nr <- dim(match$matching)[1]
    n1 <- match$matching[nr, 1]
    n2 <- match$matching[nr, 2]
    if (match$edge.length[n1] == 0) {
        match$edge.length[n1] <- match$edge.length[n2]
        match$edge.length[n2] <-  0
    }
    match$matching[nr, 3] <- 0

    if (Y$npara) para <- rep(1, Y$npara)
    else Y$pim.para <- para <- 0

    if (Y$nalpha) alpha <- rep(.5, Y$nalpha)
    else {
        Y$ncat <- rep(1, npart)
        alpha <- Y$pim.alpha <- 0
    }

    if (Y$ninvar) invar <- rep(0.5, Y$ninvar)
    else invar <- Y$pim.invar <- 0

    loglik <- 0

    if (!quiet) cat("Fitting in progress... ")
    ans <- .C("mlphylo_DNAmodel", as.integer(nb.tip), as.integer(S),
              as.double(Y$XA), as.double(Y$XC), as.double(Y$XG), as.double(Y$XT),
              as.double(Y$w), as.integer(match$matching[, 1]),
              as.integer(match$matching[, 2]), as.integer(match$matching[, 3]),
              as.double(match$edge.length), as.integer(npart),
              as.integer(Y$partition), as.integer(Y$submo), as.double(Y$xi),
              as.double(para), as.integer(Y$npara), as.integer(Y$pim.para),
              as.double(alpha), as.integer(Y$nalpha),
              as.integer(Y$pim.alpha), as.integer(Y$ncat),
              as.double(invar), as.integer(Y$ninvar),
              as.integer(Y$pim.invar), as.double(BF), as.integer(search.tree),
              as.double(loglik), NAOK = TRUE, PACKAGE = "ape")
    if (!quiet) cat("DONE!\n")
    tree <- list(matching = cbind(ans[[8]], ans[[9]], ans[[10]]),
                 edge.length = ans[[11]], tip.label = phy$tip.label)
    class(tree) <- "matching"
    tree$matching[length(tree$matching)] <-
      tree$matching[length(tree$matching) - 1] + 1
    tree <- unroot(as.phylo(tree))
    attr(tree, "loglik") <- ans[[28]]
    para <- list()
    length(para) <- npart
    names(para) <- paste("partition", 1:npart)
    for (i in 1:npart) {
        l <- list(.subst.model[ans[[14]]][i])
        names(l) <- "substitution model"
        if (Y$npara && sum(Y$pim.para[, i])) {
            tmp <- list(ans[[16]][as.logical(Y$pim.para[, i])])
            names(tmp) <- "substitution rates"
            l <- c(l, tmp)
        }
        if (Y$nalpha && sum(Y$pim.alpha[, i])) {
            tmp <- list(ans[[19]][as.logical(Y$pim.alpha[, i])])
            names(tmp) <- "alpha"
            l <- c(l, tmp)
        }
        if (Y$ninvar && sum(Y$pim.invar[, i])) {
            tmp <- list(ans[[23]][as.logical(Y$pim.invar[, i])])
            names(tmp) <- "invariants rate"
            l <- c(l, tmp)
        }
        para[[i]] <- l
    }
    attr(tree, "para") <- para
    if (dim(Y$partition)[2] > 1) attr(tree, "xi") <- ans[[15]]
    tree
}

DNAmodel <- function(model = "K80", part.model = 1,
                     ncat = 1, part.gamma = 1,
                     invar = FALSE, part.invar = 1)
{
    obj <- list(model = model, part.model = part.model,
                ncat = ncat, part.gamma = part.gamma,
                invar = invar, part.invar = part.invar)
    class(obj) <- "DNAmodel"
    obj
}

prepare.dna <- function(X, DNAmodel)
{
    ## First, check that all sequences in the list are of the
    ## same length. If OK, convert as a matrix.
    if (is.list(X)) {
        if (length(unique(unlist(lapply(X, length)))) > 1)
            stop("sequences in list must have the same lengths")
        X <- matrix(unlist(X), nrow = length(X), byrow = TRUE)
    }
    if (is.data.frame(X)) X <- as.matrix(X)
    L <- dim(X)[2]

    ## In the following 'partition' is an indicator of a set
    ## of sites that are assumed to evolve under the same
    ## parameters
    x <- rep(DNAmodel$part.model, length.out = L)
    y <- rep(DNAmodel$part.gamma, length.out = L)
    z <- rep(DNAmodel$part.invar, length.out = L)
    partition <- factor(paste(x, y, z, sep = "-"))
    npart <- nlevels(partition)

    ## find which substitution model for each partition:
    submo <- integer(npart)
    pim.submo <- matrix(0, length(DNAmodel$part.model), npart)
    for (i in 1:npart) {
        j <- as.numeric(unlist(strsplit(levels(partition)[i], "-"))[1])
        submo[i] <- which(.subst.model == DNAmodel$model[j])
        pim.submo[j, i] <- 1
    }
    ## determine how many free substitution parameters
    free.para <- c(0, 1, 0, 1, 1, 1, 2, 5)
    npara <- sum(free.para[submo])


    pim.para <- matrix(0, npara, npart)
    i <- 1
    for (j in 1:npart) {
        if (!free.para[submo[i]]) next
        pim.para[i:(i - 1 + free.para[submo[i]]), j] <- 1
        i <- i + free.para[submo[i]]
    }

    ## inter-sites variation for each partition:
    nalpha <- sum(DNAmodel$ncat > 1)
    pim.alpha <- matrix(0, nalpha, npart)
    ncat <- rep(1, npart)
    if (nalpha) {
        for (i in 1:npart) {
            j <- as.numeric(unlist(strsplit(levels(partition)[i], "-"))[2])
            pim.alpha[j, i] <- 1
            ncat[i] <- DNAmodel$ncat[j]
        }
    }

    ## proportion of invariants for each partition:
    ninvar <- sum(DNAmodel$invar)
    pim.invar <- matrix(0, ninvar, npart)
    if (ninvar) {
        for (i in 1:npart) {
            j <- as.numeric(unlist(strsplit(levels(partition)[i], "-"))[3])
            pim.invar[j, i] <- 1
        }
    }

    XA <- XC <- XG <- XT <- matrix(NA, dim(X)[1], 0)
    weight <- numeric(0)
    part <- numeric(0) # gives the length of each partition
    ## For each partition...
    for (i in 1:npart) {
        M <- X[, partition == levels(partition)[i], drop = FALSE]
        ## Get the (logical) indices of the variable sites:
        polymorph <- apply(M, 2, function(x) length(unique(x)) > 1)
        ## Frequencies of the monomorph sites:
        if (!all(polymorph)) { # is there at least one monomorph site?
            TAB <- table(M[1, !polymorph, drop = FALSE])
            x <- matrix(names(TAB), dim(M)[1], length(TAB), byrow = TRUE)
        }
        ## Frequencies of the patterns among species for the
        ## variables sites:
        if (any(polymorph)) { # is there at least one polymorph site?
            TAB2 <- table(apply(M[, polymorph, drop = FALSE], 2,
                                paste, collapse = ""))
            ## The names of TAB2 (= patterns among species) are already
            ## sorted in alphabetical order.
            y <- matrix(unlist(strsplit(names(TAB2), NULL)),
                        dim(M)[1], length(TAB2))
        }
        if (all(polymorph)) {
            dna <- y
            weight <- c(weight, as.vector(TAB2))
        } else {
            if (!any(polymorph)) {
                dna <- x
                weight <- c(weight, as.vector(TAB))
            } else {
                dna <- cbind(x, y)
                weight <- c(weight, as.vector(TAB), as.vector(TAB2))
            }
        }
        ## Transform the nucleotide values in "likelihood" values:
        xa <- xc <- xg <- xt <- matrix(0, nrow(dna), ncol(dna))
        xa[dna == "a"] <- 1
        xc[dna == "c"] <- 1
        xg[dna == "g"] <- 1
        xt[dna == "t"] <- 1
### For the moment gaps ("-") are treated in
### the same way than missing data ("n").
        uk <- dna %in% c("n", "-") # unknown
        xa[uk] <- xc[uk] <- xg[uk] <- xt[uk] <- 1
        xa[dna == "m"] <- xc[dna == "m"] <- 1
        xa[dna == "r"] <- xg[dna == "r"] <- 1
        xa[dna == "w"] <- xt[dna == "w"] <- 1
        xc[dna == "s"] <- xg[dna == "s"] <- 1
        xc[dna == "y"] <- xt[dna == "y"] <- 1
        xg[dna == "k"] <- xt[dna == "k"] <- 1
        xa[dna == "v"] <- xc[dna == "v"] <- xg[dna == "v"] <- 1
        xa[dna == "h"] <- xc[dna == "h"] <- xt[dna == "h"] <- 1
        xa[dna == "d"] <- xg[dna == "d"] <- xt[dna == "d"] <- 1
        xc[dna == "b"] <- xg[dna == "b"] <- xt[dna == "b"] <- 1

        XA <- cbind(XA, xa)
        XC <- cbind(XC, xc)
        XG <- cbind(XG, xg)
        XT <- cbind(XT, xt)
        part <- c(part, dim(xa)[2])
    }
    ## Here we add rows to these matrices for the nodes:
    tmp <- matrix(1, nrow(XA) - 2, ncol(XA))
    XA <- rbind(XA, tmp)
    XC <- rbind(XC, tmp)
    XG <- rbind(XG, tmp)
    XT <- rbind(XT, tmp)
    ## 'partition' gives the start and end of each partition:
    partition <- matrix(1, 2, npart)
    partition[2, ] <- cumsum(part)
    if (npart > 1) {
        partition[1, 2:npart] <- partition[2, 1:(npart - 1)] + 1
        partition[2, npart] <- length(weight)
        xi <- rep(1, npart - 1)
    } else xi <- 0
    list(XA = XA, XC = XC, XG = XG, XT = XT, weight = weight,
         partition = partition, submo = submo, xi = xi,
         npara = npara, pim.para = pim.para, nalpha = nalpha,
         ncat = ncat, pim.alpha = pim.alpha, ninvar = ninvar,
         pim.invar = pim.invar)
}
