### scales.R  (2004-12-18)
###
###     Add a Scale Bar or Axis to a Phylogeny Plot
###
### add.scale.bar: add a scale bar to a phylogeny plot
### axisPhylo: add a scale axis on the side of a phylogeny plot
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

add.scale.bar <- function(x = 0, y = 1, length = NULL, ...)
{
    if (is.null(length)) {
        nb.digit <- ceiling(log10(mean(.last_plot.phylo$xx))) - 2
        length <- eval(parse(text = paste("1e", nb.digit, sep = "")))
    }
    segments(x, y, x + length, y)
    text(x + length * 1.1, y, as.character(length), adj = c(0, 0.5), ...)
}

axisPhylo <- function(side = 1, ...)
{
    if (.last_plot.phylo$type %in% c("phylogram", "cladogram")) {
        if (.last_plot.phylo$direction %in% c("rightwards", "leftwards")) {
            x <- pretty(.last_plot.phylo$xx)
            if (.last_plot.phylo$direction == "rightwards")
              maxi <- max(.last_plot.phylo$xx)
            else {
                maxi <- min(.last_plot.phylo$xx)
                x <- -x
            }
        } else {
            x <- pretty(.last_plot.phylo$yy)
            if (.last_plot.phylo$direction == "upwards")
            maxi <- max(.last_plot.phylo$yy)
            else {
                maxi <- min(.last_plot.phylo$yy)
                x <- -x
            }
        }
    }
    axis(side = side, at = c(maxi - x), labels = abs(x), ...)
}
