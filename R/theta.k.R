### theta.k.R  (2002-08-28)
###
###     Population Parameter THETA using Expected Number of Alleles
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

theta.k <- function(x, n = NULL, k = NULL)
{
    if (is.null(n)) {
        if (!is.factor(x)) {
            if (is.numeric(x)) {
                n <- sum(x)
                k <- length(x)
            }
            else x <- factor(x)
        }
        if (is.factor(x)) { # ne pas remplacer par `else'...
            n <- length(x)
            k <- nlevels(x)
        }
    }
    f <- function(th) th * sum(1 / (th + (0:(n - 1)))) - k
    th <- uniroot(f, interval = c(1e-8, 100))$root
    return(th)
}
