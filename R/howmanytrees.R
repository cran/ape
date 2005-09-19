### howmanytrees.R  (2004-12-23)
###
###     Calculate Numbers of Phylogenetic Trees
###
### Copyright 2004 Emmanuel Paradis
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

howmanytrees <- function(n, rooted = TRUE, binary = TRUE,
                         labeled = TRUE, detail = FALSE)
{
    if (!labeled && !(rooted & binary))
      stop("can compute number of unlabeled trees only for rooted binary cases.")
    if (n < 3) N <- 1 else {
        if (labeled) {
            if (!rooted) n <- n - 1
            if (binary) {
                N <- prod(seq(1, (2 * n - 3), by = 2))
            } else {
                N <- matrix(0, n, n - 1)
                N[1:n, 1] <- 1
                for (i in 3:n)
                  for (j in 2:(i - 1))
                    N[i, j] <- (i + j - 2) * N[i - 1, j - 1] + j * N[i - 1, j]
                if (detail) {
                    rownames(N) <- 1:n
                    colnames(N) <- 1:(n - 1)
                } else N <- sum(N[n, ])
            }
        } else {
            N <- numeric(n)
            N[1] <- 1
            for (i in 2:n)
              if (i %% 2) N[i] <-sum(N[1:((i - 1)/2)] * N[(i - 1):((i + 1)/2)]) else {
                  x <- N[1:(i/2)]
                  y <- N[(i - 1):(i/2)]
                  y[length(y)] <- (y[length(y)] + 1) / 2
                  N[i] <- sum(x * y)
              }
            if (detail) names(N) <- 1:n else N <- N[n]
        }
    }
    N
}
