### heterozygosity.R (2002-08-28)
###
###     Heterozygosity at a Locus Using Gene Frequencies
###
### Copyright 2002 Emmanuel Paradis
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

heterozygosity <- function(x, variance = FALSE)
{
    if (!is.factor(x)) {
        if (is.numeric(x)) {
            n <- sum(x)
            k <- length(x)
            freq <- x/n
        }
        else x <- factor(x)
    }
    if (is.factor(x)) { # ne pas remplacer par `else'...
        n <- length(x)
        k <- nlevels(x)
        freq <- table(x)/n
    }
    sp2 <- sum(freq^2)
    H <- n * (1 - sp2) / (n - 1)
    if (variance) {
        sp3 <- sum(freq^3)
        var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
        return(c(H, var.H))
    }
    else return(H)
}
H <- heterozygosity
