### nodelabels.R (2006-05-31)
###
###        Labelling the Nodes and the Tips of a Tree
###
### Copyright 2004-2006 Emmanuel Paradis, 2006 Ben Bolker, and 2006 Jim Lemon
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

### from JL:
## floating.pie() from plotrix with two changes: (1) aspect ratio fixed, so pies
## will appear circular (radius is the radius in user coordinates along the x axis);
## (2) zero values allowed (but not negative)
floating.pie.asp <- function (xpos, ypos, x, edges = 200, radius = 1, col = NULL,
    startpos = 0, ...)
{
  u <- par("usr")
  user.asp <- diff(u[3:4])/diff(u[1:2])
  p <- par("pin")
  inches.asp <- p[2]/p[1]
  asp <- user.asp/inches.asp
  if (!is.numeric(x) || any(is.na(x) | x < 0)) {
    ## browser()
    stop("floating.pie: x values must be non-negative")
  }
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  if (is.null(col))
    col <- rainbow(nx)
  else if (length(col) < nx)
    col <- rep(col, nx)
  bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
  for (i in 1:nx) {
    n <- max(2, floor(edges * dx[i]))
    t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
    xc <- c(cos(t2p) * radius + xpos, xpos)
    yc <- c(sin(t2p) * radius*asp + ypos, ypos)
    polygon(xc, yc, col = col[i], ...)
    ## t2p <- 2 * pi * mean(x[i + 0:1]) + startpos
    ## xc <- cos(t2p) * radius
    ## yc <- sin(t2p) * radius*asp
    ## lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
  }
  ## return(bc)
}

BOTHlabels <- function(text, sel, adj, frame, pch, thermo,
                       pie, piecol, col, bg, ...)
{
    if (missing(text)) text <- NULL
    if (length(adj) == 1) adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie))
      text <- sel
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
        parusr <- par("usr")
        width <- CEX * (parusr[2] - parusr[1]) / 40
        height <- CEX * (parusr[4] - parusr[3]) / 15
        if (is.vector(thermo)) thermo <- cbind(thermo, 1 - thermo)
        thermo <- height * thermo
        xl <- .last_plot.phylo$xx[sel] - width/2
        xr <- xl + width
        yc <- .last_plot.phylo$yy[sel]
        yb <- yc - height/2
        if (is.null(piecol)) piecol <- rainbow(ncol(thermo))
        ## draw the first rectangle:
        rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
        for (i in 2:ncol(thermo))
          rect(xl, yb + rowSums(thermo[, 1:(i - 1), drop = FALSE]),
               xr, yb + rowSums(thermo[, 1:i]),
               border = NA, col = piecol[i])
        rect(xl, yb, xr, yb + height, border = "black")
        segments(xl, yc, xl - width/5, yc)
        segments(xr, yc, xr + width/5, yc)
    }
    ## from BB:
    if (!is.null(pie)) {
        if (is.vector(pie)) pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2]) / 50
        for (i in 1:length(sel))
          floating.pie.asp(.last_plot.phylo$xx[sel[i]],
                           .last_plot.phylo$yy[sel[i]],
                           pie[i, ], radius = xrad, col = piecol)
    }
    if (!is.null(text)) text(.last_plot.phylo$xx[sel],
                             .last_plot.phylo$yy[sel],
                             text, adj = adj, col = col, ...)
    if (!is.null(pch)) points(.last_plot.phylo$xx[sel] + adj[1] - 0.5,
                              .last_plot.phylo$yy[sel] + adj[2] - 0.5,
                              pch = pch, col = col, bg = bg, ...)
}

nodelabels <- function(text, node, adj = c(0.5, 0.5), frame = "rect",
                       pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
                       col = "black", bg = "lightblue", ...)
{
    sel <- if (missing(node))
      names(.last_plot.phylo$xx)[as.numeric(names(.last_plot.phylo$xx)) < 0]
    else as.character(-abs(as.numeric(node)))
    BOTHlabels(text, sel, adj, frame, pch, thermo, pie, piecol, col, bg, ...)
}

tiplabels <- function(text, tip, adj = c(0.5, 0.5), frame = "rect",
                      pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
                      col = "black", bg = "yellow", ...)
{
    sel <- if (missing(tip))
      names(.last_plot.phylo$xx)[as.numeric(names(.last_plot.phylo$xx)) > 0]
    else tip
    BOTHlabels(text, sel, adj, frame, pch, thermo, pie, piecol, col, bg, ...)
}
