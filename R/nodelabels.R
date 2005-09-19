### nodelabels.R (2005-09-12)
###
###             Labelling the Nodes of a Tree
###
### Copyright 2004-2005 Emmanuel Paradis
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

nodelabels <- function(text, node, adj = c(0.5, 0.5), frame = "rect",
                       pch = NULL, thermo = NULL, col = "black",
                       bg = "lightblue", ...)
{
    sel <- if (missing(node))
      names(.last_plot.phylo$xx)[as.numeric(names(.last_plot.phylo$xx)) < 0]
    else as.character(-abs(as.numeric(node)))
    if (missing(text)) text <- NULL
    if (length(adj) == 1) adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo))
      text <-
        names(.last_plot.phylo$xx)[as.numeric(names(.last_plot.phylo$xx)) < 0]
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) args$cex else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            offs <- xinch(0.03)
            xl <- .last_plot.phylo$xx[sel] - strwidth(text) * CEX * adj[1] - offs
            xr <- xl + strwidth(text) * CEX + 2 * offs
            yb <- .last_plot.phylo$yy[sel] - strheight(text) * CEX * adj[2] - offs
            yt <- yb + strheight(text) * CEX + 2 * offs
            rect(xl, yb, xr, yt, col = bg)
        }
        if (frame == "circle") {
            radii <- apply(cbind(strheight(text), strwidth(text)), 1, max) * 0.6
            symbols(.last_plot.phylo$xx[sel], .last_plot.phylo$yy[sel],
                    circles = radii, inches = FALSE, add = TRUE,
                    bg = bg)
        }
    }
    if (!is.null(thermo)) {
        width <- CEX * (par("usr")[2] - par("usr")[1]) / 30
        height <- CEX * (par("usr")[4] - par("usr")[3]) / 15
        symbols(.last_plot.phylo$xx[sel], .last_plot.phylo$yy[sel],
                thermometers = cbind(width, height, thermo),
                inches = FALSE, add = TRUE, fg = col, bg = bg)
    }
    if (!is.null(text)) text(.last_plot.phylo$xx[sel],
                             .last_plot.phylo$yy[sel],
                             text, adj = adj, ...)
    if (!is.null(pch)) points(.last_plot.phylo$xx[sel],
                              .last_plot.phylo$yy[sel],
                              pch = pch, col = col, bg = bg, ...)
}
