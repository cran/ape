### dist.dna.R (2005-09-12)
###
###     Pairwise Distances from DNA Sequences
###
### Copyright 2002-2005 Emmanuel Paradis
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

dist.dna <- function(x, model = "K80", variance = FALSE,
                     gamma = FALSE, pairwise.deletion = FALSE,
                     base.freq = NULL, as.matrix = FALSE)
{
    MODELS <- c("raw", "JC69", "K80", "F81", "K81",
                "F84", "T92", "TN93", "GG95")
    imod <- which(MODELS == model)
    if (gamma && imod %in% c(3, 5, 6, 8)) {
        warning(paste("gamma-correction not available for model", model))
        gamma <- FALSE
    }
    ## if the data are in matrix-form, transpose them to have
    ## individuals as cols and sites as rows
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) x <- t(x)
    ## if the data are given as a list, we transform it as a matrix
    ## with individuals as cols and sites as rows
    if (is.list(x)) {
        nm <- names(x)
        n <- length(x)
        if (length(unique(unlist(lapply(x, length)))) != 1)
            stop("DNA sequences in list not of the same length.")
        x <- unlist(x)
        nL <- length(x)
        dim(x) <- c(nL / n, n)
        colnames(x) <- nm
    }
    BF <- if (is.null(base.freq)) base.freq(x) else base.freq
    if (!pairwise.deletion) {
        sel <- !apply(x, 1, function(x) any(x == "n"))
        x <- x[sel, ]
    }
    s <- nrow(x) # the number of sites
    n <- ncol(x) # the number of individuals
    var <- if (variance) numeric(n * (n - 1) / 2) else 0
    if (!gamma) gamma <- alpha <- 0
    else {
        alpha <- gamma
        gamma <- 1
    }
    d <- .C("dist_dna", as.character(x), as.integer(n), as.integer(s),
            as.integer(imod), as.double(numeric(n * (n - 1) / 2)),
            as.double(BF), as.integer(pairwise.deletion),
            as.integer(variance), as.double(var),
            as.integer(gamma), as.double(alpha),
            NAOK = TRUE, PACKAGE = "ape")
    if (variance) var <- d[[9]]
    d <- d[[5]]
    attr(d, "Size") <- n
    attr(d, "Labels") <- colnames(x)
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- model
    if (variance) attr(d, "variance") <- var
    class(d) <- "dist"
    if (as.matrix) d <- as.matrix(d)
    d
}
