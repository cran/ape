### seg.sites.R  (2004-03-19)
###
###     Find Segregating Sites in DNA Sequences
###
### Copyright 2004 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

seg.sites <- function(X)
{
    if (is.list(X)) {
        if (length(unique(unlist(lapply(X, length)))) > 1)
          stop("sequences in list must have the same lengths")
        X <- matrix(unlist(X), nrow = length(X), byrow = TRUE)
    }
    if (is.data.frame(X)) X <- as.matrix(X)
    which(apply(X, 2, function(x) length(unique(x)) > 1))
}
