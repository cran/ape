### base.freq.R  (2002-08-28)
###
###     Base frequencies from DNA Sequences
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

base.freq <- function(x)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) dim(x) <- NULL
    if (is.list(x)) x <- unlist(x)
    x <- x[grep("[acgt]", x)] # get only the known bases
    n <- length(x)
    return(table(x) / n)
}
