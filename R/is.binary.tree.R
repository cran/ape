### is.binary.tree.R  (2002-09-12)
###
###     Tests whether a given phylogenetic tree is binary
###
### Copyright 2002 Korbinian Strimmer <strimmer@stat.uni-muenchen.de>
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

is.binary.tree <- function(tree)
{
    if (class(tree) != "phylo") stop("object \"tree\" is not of class \"phylo\"")
    ## modified by EP so that it works without edge lengths too (2005-05-31):
    tmp <- as.numeric(tree$edge)
    if (max(tmp) - 1 ==  -min(tmp)) return(TRUE)
    else return(FALSE)
}
