### zoom.R  (2003-12-27)
###
###     Zoom on a Portion of a Phylogeny
###
### Copyright 2003 Emmanuel Paradis <paradis@isem.univ-montp2.fr>
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

zoom <- function(phy, focus, type = "phylogram", use.edge.length = TRUE)
{
    if (is.character(focus)) focus <- which(phy$tip.label == focus)
    ext <- drop.tip(phy, phy$tip.label[-focus])
    layout(matrix(1:2, 1, 2), c(1, 3))
    phy$tip.label[-focus] <- ""
    phy$tip.label[focus] <- "o"
    par(col = "blue")
    plot.phylo(phy, type = type, use.edge.length = use.edge.length,
               label.offset = 0, no.margin = TRUE)
    par(col = "black")
    plot.phylo(ext, type = type, use.edge.length = use.edge.length,
               label.offset = 0, no.margin = TRUE)
}
