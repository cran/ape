### cherry.R  (2002-08-28)
###
###     Number of Cherries and Null Models of Trees
###
### Copyright 2002 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

cherry <- function(phy)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    tmp <- as.numeric(phy$edge)
    n <- max(tmp)
    nb.node <- -min(tmp)
    if (nb.node != n - 1) stop("\"phy\" is not fully dichotomous")
    if (n < 4) stop("not enough tips in your phylogeny for this analysis")	
    cherry <- sum(table(phy$edge[, 1][phy$edge[, 2] > 0]) == 2)
    if (n < 20) small.n <- TRUE else small.n <- FALSE
    if (small.n) {
        P.yule <- f.cherry.yule(n, cherry)
        P.uniform <- f.cherry.uniform(n, cherry)
    }
    else {
        P.yule <- 2 * (1 - pnorm(abs(cherry - n / 3) / sqrt(2 * n / 45)))
        mu.unif <- n * (n - 1) / (2 * (2 * n - 5))
        sigma2.unif <- n * (n - 1) * (n - 4) * (n - 5) / (2 * (2 * n - 5)^2 * (2 * n -7))
        P.uniform <- 2 * (1 - pnorm(abs(cherry - mu.unif) / sqrt(sigma2.unif)))
    }
    cat("\nAnalysis of the Number of Cherries in a Tree\n\n")
    cat("Phylogenetic tree:", deparse(substitute(phy)), "\n")
    cat("Number of tips:", n, "\n")
    cat("Number of cherries:", cherry, "\n\n")
    cat("Null hypothesis: Yule model\n")
    cat("    P-value =", round(P.yule, 4), "\n\n")
    cat("Null hypothesis: uniform model\n")
    cat("    P-value =", round(P.uniform, 4), "\n\n")
    if (!small.n) cat("(P-values were computed using normal approximations)\n")
}

f.cherry.yule <- function(n, k)
{
    if (k == 0 | k > floor(n/2)) P <- 0
    else {
        if (n == 4) {
            if (k == 1) P <- 2/3 else if (k == 2) P <- 1/3 else P <- 0
        }
        else {
            P <- (1 - 2 * (k - 1)/(n - 1)) * f.cherry.yule(n - 1, k - 1) +
              2 * k/(n - 1) * f.cherry.yule(n - 1, k)
        }
    }
    return(P)
}

f.cherry.uniform <- function(n, k)
{
    if (k == 0 | k > floor(n/2)) P <- 0
    else {
        if (n == 4) {
            if (k == 1) P <- 4/5 else if (k == 2) P <- 1/5 else P <- 0
        }
        else {
            if (k == 1) P <- 0
            else
              P <- (gamma(n + 1) * gamma(n - 2 + 1) * gamma(n - 4 + 1) * 2^(n-2*k)) /
                (gamma(n -2 * k + 1) * gamma(2 * n - 4 + 1) * gamma(k + 1)
                 * gamma(k - 2 + 1))
        }
    }
    return(P)
}
