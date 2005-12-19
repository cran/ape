### MoranI.R  (2004-09-16)
###
###     Moran's I Autocorrelation Index
###
### Copyright 2004 Julien Dutheil <julien.dutheil@univ-montp2.fr>
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

Moran.I <- function(
  x,              # x the vector of values to analyse
  dist,           # dist is the distance matrix to use
  scaled = FALSE, # tell if we must scale the indice to allow comparisons
  na.rm = FALSE   # should we ignore missing values?
) {
  # Number of values:
  nn <- length(x);
  n <- ifelse(na.rm, sum(!is.na(x)), nn)
  if(dim(dist)[1] != nn | dim(dist)[2] != nn) {
    stop("\"dist\" must be a matrix of size n*n, with n=length(x)=", n, ".");
  }

  # Normalization of distances:
  for(i in 1:nn) {
    s <- sum(dist[i,]);
    if(s != 0) {
      dist[i,] <- dist[i,] / s;
    } # Else row is only made of zeros, need to use scale=T for comparisons.
  }

  s <- sum(dist);
  m <- mean(x, na.rm = na.rm);
  c <- 0; # covariance
  v <- 0; # variance
  # Compute covariance:
  for(i in 1:nn) {
    if(na.rm && is.na(x[i])) next
    for(j in 1:nn) {
      if(na.rm && is.na(x[j])) next
      c <- c + dist[i,j] * (x[i] - m) * (x[j] - m);
    }
  }
  # Compute variance:
  for(i in 1:nn) {
    if(na.rm && is.na(x[i])) next
    v <- v + (x[i] - m) ^ 2;
  }
  obs <- (n / s) * (c / v);

  # Now computes expected mean and sd
  # under the randomization null hypothesis:

  # Expected mean:
  ei <- -1/(n - 1);

  # Scaling:
  if(scaled) {
    temp <- vector(length = nn);
    for(i in 1:nn) {
      temp[i] <- (x[i] - m) * sum(dist[i,]);
    }
    i.max <- (n/s) * (sd(temp, na.rm = na.rm) / sd(x - m, na.rm = na.rm));
    obs <- obs/i.max;
  }

  # Expected sd:
  S1 <- (1/2) * sum((dist + t(dist))^2);
  S2 <- 0;
  for(i in 1:nn) {
    S2 <- S2 + (sum(dist[i,]) + sum(dist[,i]))^2;
  }
  k <- ((1/n) * sum((x - m)^4, na.rm = na.rm)) / ((1/n) * sum((x - m)^2, na.rm = na.rm))^2
  sdi <- sqrt((
          n * ((n^2 - 3*n + 3)*S1 - n*S2 + 3*s^2)
          - k * (n*(n-1)*S1 - 2*n*S2 + 6*s^2)
         ) /
         ((n-1)*(n-2)*(n-3)*s^2) -
         1/((n-1)^2));

  names(obs) <- NULL
  # Computes p-value:
  pv <- 1 - 2*abs(pnorm(obs, m = ei, sd = sdi) - 0.5);

  return(list(observed = obs, expected = ei, sd = sdi, p.value = pv));
}

weight.taxo <- function (x) {
  n <- length(x)
  d <- matrix(ncol = n, nrow = n, 0)
  for (i in 1:n) {
    d[i, which(x[i] == x)] <- 1
    d[i, i] <- 0
  }
  return(d)
}

correlogram.formula <- function(formula, data = NULL, use = "all.obs")
{
  err <- "Formula must be of the kind \"y1+y2+..+yn~x1/x2/../xn\"."
  if(!(use %in% c("all.obs", "complete.obs", "pairwise.complete.obs")))
    stop("Argument 'use' must be either 'all.obs', 'complete.obs' or 'pairwise.complete.obs'.")

  if (is.null(data)) data <- parent.frame()

  if(formula[[1]] != "~") stop(err);

  #Must check if y is transformed:
  get.var <- function(var) {
    if(length(var) == 1) {
      # Simple variable
      var.name <- deparse(var)
      if(!is.null(data[[var.name]])) {
        # Look within dataframe:
		    y <- data[[var.name]]
      } else {
        # Not found, look in global environment:
		    y <- get(var.name)
      }
      return(list(y = y, name = var.name))
    } else if(length(var) == 2) {
      # Transformed variable:
      var.name <- deparse(var[[2]])
      fun.name <- deparse(var[[1]])
      if(!is.null(data[[var.name]])) {
        # Look within dataframe:
		    y <- data[[var.name]]
      } else {
        # Not found, look in global environment:
		    y <- parent.frame(2)[[var.name]]
      }
      return(list(y=get(fun.name)(y), name=deparse(var)))
    } else if (length(var) == 3) {
      if(var[[1]] == '$') {
        var.name <- deparse(var[[3]])
        df.name  <- deparse(var[[2]])
        return(list(y = get(df.name)[[var.name]], name = deparse(var)))
      } else stop(err)
    }
  }

  y <- list()
  ally <- formula[[2]]
  while(length(ally) == 3 && ally[[1]] == '+') {
    var <- get.var(ally[[3]])
    y[[var$name]] <- var$y
    ally <- ally[[2]]
  }
  # Last y:
  var <- get.var(ally)
  y[[var$name]] <- var$y

  ##Groups:
  groups <- formula[[3]]
  d <- list()
  g <- list()

  while(length(groups) == 3 && groups[[1]] == '/') {
    group  <- get.var(groups[[3]])
    groups <- groups[[2]]
    cat("Analysing level:", group$name, "\n")
    g[[group$name]] <- group$y
    d[[group$name]] <- weight.taxo(group$y)
  }
  # The last group:
  group <- get.var(groups)
  #cat("Analysing level:", group$name, "\n")
  g[[group$name]] <- group$y
  d[[group$name]] <- weight.taxo(group$y)

  # Remove all data with missing grouping values:
  filter <- rep(TRUE, length(y[[1]])) # All obs used
  if(use == "complete.obs" || use == "pairwise.complete.obs") {
    G <- sapply(g, is.na)
    for(i in 1:dim(G)[[2]]) {
      filter <- filter & !G[,i]
    }
  }

  # Deal with the complete.obs option:
  if(use == "complete.obs") {
    M <- sapply(y, is.na)
    for(i in 1:dim(M)[[2]]) {
      filter <- filter & !M[,i]
    }
  }

  # Now compute Moran's I:
  n <- length(d)
  l <- p <- i <- vector(length = n)
  corList <- list()
  for(k in names(y)) {
    for(j in 1:n) {
      if(j == 1) Mat <- d[[j]]
      else Mat <- d[[j]] & !d[[j-1]]
      I.M  <- Moran.I(y[[k]][filter], Mat[filter, filter], scale = TRUE, na.rm = (use == "pairwise.complete.obs"));
      i[j] <- I.M$obs
      p[j] <- I.M$p.v
      l[j] <- names(d)[j]
    }

    # Create an object of class 'correlogram':
    corr <- list(obs = i, p.values = p, labels = l, filter = filter)
    class(corr) <- "correlogram"
    corList[[k]] <- corr
  }
  class(corList) <- "correlogramList"
  if(length(corList) == 1) return(corList[[names(y)]])
  else                     return(corList)
}

discrete.dist <- function(dist, inf, sup)
{
  if (class(dist) != "matrix") stop("object \"dist\" is not of class \"matrix\"")
  n <- dim(dist)[1]
  d <- matrix(ncol = n, nrow = n)
  rownames(d) <- rownames(dist)
  colnames(d) <- colnames(dist)
  for(i in 1:n) {
    for(j in 1:n) {
      d[i, j] <- ifelse(dist[i,j] > inf & dist[i,j] <= sup & i != j, 1, 0)
    }
  }
  return(d)
}

correlogram.phylo <- function(x, phy, nclass = NULL, breaks = NULL)
{
  if (!("phylo" %in% class(phy))) stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) stop("tree \" phy\" must have branch lengths.")
  #Get the minimum and maximum distance in the tree:
  dist <- cophenetic.phylo(phy)
  #What classes to use?
  if(!is.null(breaks)) {
    # User-defined breaks:
    s <- sort(breaks)
    nclass <- length(s)-1
  } else if(!is.null(nclass)) {
    # Equal classes:
    l.min <- min(dist[dist != 0])
    l.max <- max(dist)
    if(nclass < 2) stop("\"nclass\" must be > 1.")
    s <- seq(from=l.min, to=l.max, length=nclass+1)
  }
  if(!is.null(nclass)) {
    l <- p <- i <- vector(length = nclass)
    for(j in 1:(nclass)) {
      cat("Analysing level:",j,"\n");
      Mat <- discrete.dist(dist, s[j], s[j+1])
      I.M <- Moran.I(x, Mat, scale=TRUE)
      i[j] <- I.M$obs
      p[j] <- I.M$p.v
      l[j] <- paste("Class",j)
    }
  } else {
    cat("using whole matrix\n");
    I.M <- Moran.I(x, dist)
    i <- I.M$obs
    p <- I.M$p.v
    l <- paste("Distance")
  }

  # Create an object of class 'correlogram':
  corr <- list(obs=i, p.values=p, labels=l)
  class(corr) <- "correlogram"
  return(corr)
}

plot.correlogram <- function(x, test.level=0.05, ...)
{
  if (!("correlogram" %in% class(x))) stop("object \"x\" is not of class \"correlogram\"")
  # Black circles are significant at the 5% level:
  pch <- ifelse(x$p.values < test.level, 19, 21)
  # Plot it!
  return(xyplot(x$obs~ordered(x$l,levels=x$l), type="b", xlab="Rank", ylab="I / Imax", lty=2, lwd=2, cex=1.5, pch=pch, ...))
}

panel.superpose.correlogram <- function (x, y = NULL, subscripts, groups, panel.groups = "panel.xyplot",
    col, col.line = superpose.line$col, col.symbol = superpose.symbol$col,
    pch = superpose.symbol$pch, p.values = NULL, test.level=0.05, cex = superpose.symbol$cex, font = superpose.symbol$font,
    fontface = superpose.symbol$fontface, fontfamily = superpose.symbol$fontfamily,
    lty = superpose.line$lty, lwd = superpose.line$lwd, ...)
{
    x <- as.numeric(x)
    if (!is.null(y))
        y <- as.numeric(y)
    if (length(x) > 0) {
        if (!missing(col)) {
            if (missing(col.line))
                col.line <- col
            if (missing(col.symbol))
                col.symbol <- col
        }
        superpose.symbol <- trellis.par.get("superpose.symbol")
        superpose.line <- trellis.par.get("superpose.line")
        vals <- if (is.factor(groups))
            levels(groups)
        else sort(unique(groups))
        nvals <- length(vals)
        col.line <- rep(col.line, length = nvals)
        col.symbol <- rep(col.symbol, length = nvals)
        if(is.null(p.values))
          pch <- rep(pch, length = nvals)
        else
          pch <- ifelse(p.values < test.level, 19, 21)
        lty <- rep(lty, length = nvals)
        lwd <- rep(lwd, length = nvals)
        cex <- rep(cex, length = nvals)
        font <- rep(font, length = nvals)
        fontface <- rep(fontface, length = nvals)
        fontfamily <- rep(fontfamily, length = nvals)
        panel.groups <- if (is.function(panel.groups))
            panel.groups
        else if (is.character(panel.groups))
            get(panel.groups)
        else eval(panel.groups)
        for (i in seq(along = vals)) {
            id <- (groups[subscripts] == vals[i])
            if (any(id)) {
                args <- list(x = x[id], groups = groups, subscripts = subscripts[id],
                  pch = pch[id], cex = cex[i], font = font[i],
                  fontface = fontface[i], fontfamily = fontfamily[i],
                  col.line = col.line[i], col.symbol = col.symbol[i],
                  lty = lty[i], lwd = lwd[i], ...)
                if (!is.null(y))
                  args$y <- y[id]
                do.call("panel.groups", args)
            }
        }
    }
}

plot.correlogramList <- function(x, test.level=0.05, ...)
{
  if (!("correlogramList" %in% class(x))) stop("object \"x\" is not of class \"correlogramList\"")
  #Build a dataframe:
  obs <- numeric(0)
  lev <- numeric(0)
  cor <- numeric(0)
  pvl <- numeric(0)
  for(i in names(x)) {
    obs <- c(obs, x[[i]]$obs)
    lev <- c(lev, x[[i]]$l)
    cor <- c(cor, rep(i, length(x[[i]]$obs)))
    pvl <- c(pvl, x[[i]]$p.values)
  }
  lev <- ordered(lev, levels=unique(lev))
  return(xyplot(obs~lev, groups=cor,
        type="b", xlab="Rank", ylab="I / Imax",
        lty=2, lwd=2, cex=1.5, panel=panel.superpose.correlogram,
        p.values=pvl, key=simpleKey(names(x), lines=TRUE, points=FALSE, rectangle=FALSE), ...))
}

#data(carnivora)
#co <- correlogram.formula(log10(SW) + log10(FW) ~ Order/SuperFamily/Family/Genus, data=carnivora)
#plot(co)


