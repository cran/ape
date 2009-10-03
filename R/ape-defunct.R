klastorin <- function(phy)
    .Defunct(msg = '\'klastorin\' has been removed from ape,
    see help("ape-defunct") for details.')

mlphylo <-
    function(x, phy, model = DNAmodel(), search.tree = FALSE,
             quiet = FALSE, value = NULL, fixed = FALSE)
    .Defunct(msg = '\'mlphylo\' has been removed from ape,
    see help("ape-defunct") for details.')

DNAmodel <- function(model = "K80", partition = 1,
         ncat.isv = 1, invar = FALSE,
         equal.isv = TRUE, equal.invar = 1)
    .Defunct(msg = '\'DNAmodel\' has been removed from ape,
    see help("ape-defunct") for details.')

sh.test <- function(..., x, model = DNAmodel(), B = 100)
    .Defunct(msg = '\'sh.test\' has been removed from ape,
    see help("ape-defunct") for details.')

heterozygosity <- function (x, variance = FALSE)
    .Defunct(msg = '\'heterozygosity\' has been moved from ape to pegas,
    see help("ape-defunct") for details.')

H <- function(x, variance = FALSE)
    heterozygosity (x, variance = FALSE)

nuc.div <- function(x, variance = FALSE, pairwise.deletion = FALSE)
    .Defunct(msg = '\'nuc.div\' has been moved from ape to pegas,
    see help("ape-defunct") for details.')

theta.h <- function(x, standard.error = FALSE)
    .Defunct(msg = '\'theta.h\' has been moved from ape to pegas,
    see help("ape-defunct") for details.')

theta.k <- function(x, n = NULL, k = NULL)
    .Defunct(msg = '\'theta.k\' has been moved from ape to pegas,
    see help("ape-defunct") for details.')

theta.s <- function(s, n, variance = FALSE)
    .Defunct(msg = '\'theta.s\' has been moved from ape to pegas,
    see help("ape-defunct") for details.')
