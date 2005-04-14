### nuc.div.R  (2005-04-01)
###
###     Nucleotide Diversity
###
### Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

nuc.div <- function(x, variance = FALSE, pairwise.deletion = FALSE)
{
    if (pairwise.deletion & variance)
      warning("cannot compute the variance of nucleotidic diversity
with pairwise deletion: try 'pairwise.deletion = FALSE' instead.")
    if (is.list(x)) {
        if (length(unique(unlist(lapply(x, length)))) > 1)
          stop("sequences in list must have the same lengths")
        x <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
    }
    if (is.data.frame(x)) x <- as.matrix(x)
    sum.pi <- 0
    N <- dim(x)[1] # number of sequences
    if (pairwise.deletion) {
        for (i in 1:(N - 1))
          for (j in (i + 1):N) {
              sel <- !(x[i, ] == "n" | x[j, ] == "n")
              sum.pi <- sum.pi + sum(x[i, ][sel] != x[j, ][sel]) / sum(sel)
          }
        obj <- sum.pi / (N * (N - 1) / 2)
    } else {
        sel <- !apply(x, 2, function(x) any(x == "n"))
        n <- sum(sel) # length of the sequences
        for (i in 1:(N - 1))
          for (j in (i + 1):N)
            sum.pi <- sum.pi + sum(x[i, ][sel] != x[j, ][sel]) / n
        obj <- sum.pi / (N * (N - 1) / 2)
        if (variance) {
            var <- (N + 1) * obj / (3 * (N + 1) * n) +
              2 * (N^2 + N + 3) * obj / (9 * N * (N - 1))
            obj <- c(obj, var)
        }
    }
    obj
}
