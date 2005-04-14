### ltt.plot.R  (2004-08-31)
###
###     Lineages Through Time Plot
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

ltt.plot <- function(phy, xlab = "Time", ylab = "N", ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    time <- sort(branching.times(phy))
    N <- 1:(length(time) + 1)
    plot(-c(rev(time), 0), N, xlab = xlab, ylab = ylab,
         xaxs = "r", yaxs = "r", type = "S", ...)
}

ltt.lines <- function(phy, ...)
{
    if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
    time <- sort(branching.times(phy))
    N <- 1:(length(time) + 1)
    lines(-c(rev(time), 0), N, type = "S", ...)
}

mltt.plot <- function(phy, ..., dcol = TRUE, dlty = FALSE, legend = TRUE,
                      xlab = "Time", ylab = "N")
{
    ## this will also accept objects of class `c("phylo", "multi.tree")'
    if (class(phy)[1] != "phylo")
      stop("object \"phy\" is not of class \"phylo\"")
    ltt.xy <- function(phy) {
        x <- -c(rev(sort(branching.times(phy))), 0)
        names(x) <- NULL
        y <- 1:length(x)
        cbind(x, y)
    }
    if (length(class(phy)) == 1) {
        TREES <- list(ltt.xy(phy))
        names(TREES) <- deparse(substitute(phy))
    } else {
        TREES <- lapply(phy, ltt.xy)
        names(TREES) <- names(phy)
    }
    dts <- list(...)
    if (length(dts)) {
        mc <- as.character(match.call())[-(1:2)]
        nms <- mc[1:length(dts)]
        for (i in 1:length(dts)) {
            if (length(class(dts[[i]])) == 1) {
                a <- list(ltt.xy(dts[[i]]))
                names(a) <- nms[i]
            } else {
                a <- lapply(dts[[i]], ltt.xy)
                names(a) <- names(dts[[i]])
            }
            TREES <- c(TREES, a)
        }
    }
    n <- length(TREES)
    xl <- c(min(unlist(lapply(TREES, function(x) min(x[, 1])))), 0)
    yl <- c(1, max(unlist(lapply(TREES, function(x) max(x[, 2])))))

    plot(0, 0, type = "n", xlim = xl, ylim = yl, xaxs = "r", yaxs = "r",
         xlab = xlab, ylab = ylab)

    lty <- if (!dlty) rep(1, n) else 1:n
    col <- if (!dcol) rep(1, n) else topo.colors(n)

    for (i in 1:n) lines(TREES[[i]], col = col[i], lty = lty[i], type = "S")

    if (legend) legend(xl[1], yl[2], legend = names(TREES), lty = lty, col = col, bty = "n")
}
