### add.scale.bar.R  (2002-10-02)
###
###     Add a Scale Bar to a Phylogeny Plot
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

add.scale.bar <- function(x, y, phy, bar = NULL)
{
    if (is.null(bar)) {
        nb.digit <- ceiling(log10(mean(phy$edge.length)))
        bar <- eval(parse(text = paste("1e", nb.digit - 1, sep = "")))
    }
    segments(x, y, x + bar, y)
    text(x + bar * 1.1, y, as.character(bar), adj = c(0, 0.5))
}
