### theta.s.R  (2002-08-28)
###
###     Population Parameter THETA using Segregating Sites
###       in DNA Sequences
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

theta.s <- function(s, n, variance = TRUE)
{
    a1 <- sum(1 / (1:(n - 1)))
    th <- s / a1
    if (variance) {
        a2 <- sum(1 / (1:(n - 1))^2)
        var.th <- (a1^2 * s + a2 * s^2) / (a1^2 * (a1^2 + a2))
        return(c(th, var.th))
    }
    else return(th)
}
