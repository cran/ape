### theta.h.R  (2002-08-28)
###
###     Population Parameter THETA using Homozigosity
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

theta.h <- function(x, standard.error = FALSE)
{
    HE <- H(x, variance = TRUE)
    sdH <- HE[2]
    HE <- HE[1]
    f <- function(th) HE - th * (1 + (2 * (1 + th)) / ((2 + th) * (3 + th)))
    th <- uniroot(f, interval = c(0, 1))$root
    if (standard.error) {
        SE <- (2 + th)^2 * (2 + th)^3 * sdH /
          HE^2 * (1 + th) * ((2 + th) * (3 + th) * (4 + th) + 10 * (2 + th) + 4)
        return(c(th, SE))
    }
    else return(th)
}
