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

# For debugging:
#setwd("Dev/r/APE/Data/")           
#carn<-read.table("Gittleman1986.csv", header=T, sep="\t")


Moran.I <- function(
  x,             # x the vector of values to analyse
  dist,          # dist is the distance matrix to use
  scaled = FALSE # tell if we must scale the indice to allow comparisons
) {

  # Number of values:
  n <- length(x);
  if(dim(dist)[1] != n | dim(dist)[2] != n) {
    stop("\"dist\" must be a matrix of size n*n, with n=length(x)=", n, ".");
  }

  # Normalization of distances:
  for(i in 1:n) {
    s <- sum(dist[i,]);
    if(s != 0) {
      dist[i,] <- dist[i,] / s;
    } # Else row is only made of zeros, need to use scale=T for comparisons.
  }

  s <- sum(dist);
  m <- mean(x);
  c <- 0; # covariance
  v <- 0; # variance
  # Compute covariance:
  for(i in 1:n) {
    for(j in 1:n) {
      c <- c + dist[i,j] * (x[i] - m) * (x[j] - m);
    }
  }
  # Compute variance:
  for(i in 1:n) {
    v <- v + (x[i] - m) ^ 2;
  }
  obs <- (n / s) * (c / v);

  # Now computes expected mean and sd
  # under the randomization null hypothesis:

  # Expected mean:
  ei <- -1/(n - 1);

  # Scaling:
  if(scaled) {
    temp <- vector(length=n);
    for(i in 1:n) {
      temp[i] <- (x[i] - m) * sum(dist[i,]);
    }
    i.max <- (n/s) * (sd(temp) / sd(x - m));
    obs <- obs/i.max;
  }

  # Expected sd:
  S1 <- (1/2) * sum((dist + t(dist))^2);
  S2 <- 0;
  for(i in 1:n) {
    S2 <- S2 + (sum(dist[i,]) + sum(dist[,i]))^2;
  }
  k <- ((1/n) * sum((x - m)^4)) / ((1/n) * sum((x - m)^2))^2
  sdi <- sqrt((
          n * ((n^2 - 3*n + 3)*S1 - n*S2 + 3*s^2)
          - k * (n*(n-1)*S1 - 2*n*S2 + 6*s^2)
         ) /
         ((n-1)*(n-2)*(n-3)*s^2) -
         1/((n-1)^2));

  # Computes p-value:
  pv <- 1 - 2*abs(pnorm(obs, m=ei, sd=sdi) - 0.5);

  return(list(observed = obs, expected=ei, sd=sdi, p.value=pv));
}

dist.taxo <- function(x)
{
  n <- length(x);
  d <- matrix(ncol = n, nrow = n);
  for(i in 1:n) {
    for(j in 1:n) {
      d[i,j] <- ifelse(x[i] == x[j] & i != j, 1, 0);
    }
  }
  return(d);
}

correlogram.formula <- function(formula, data)
{
  err <- "Formula must be of the kind \"y~x1/x2/../xn\"."

  if (is.null(data)) data <- parent.frame()

  if(formula[[1]] != "~") stop(err);

  # Variable:
  var <- formula[[2]]
  #Must check if y is transformed:
  if(length(var) == 1) {
    # Simple variable
    var.name <- deparse(var)
    y <- data[[var.name]]
  } else if(length(var) == 2) {
    # Transformed variable
    var.name <- deparse(var[[2]])
    fun.name <- deparse(var[[1]])
    y <- get(fun.name)(data[[var.name]])
  } else stop(err)
  #Groups:
  groups <- formula[[3]]
  d <- list()

  while(length(groups) == 3) {
    if(groups[[1]] != "/") stop(err)
    group  <- groups[[3]]
    groups <- groups[[2]]
    if(length(group) != 1) stop(err)
    s <- deparse(group)
    cat("Analysing level:", s, "\n")
    d[[s]] <- dist.taxo(data[[s]])
  }
  # The last group:
  group <- groups
  if(length(group) != 1) stop(err)
  s <- deparse(group)
  cat("Analysing level:",s,"\n")
  d[[s]] <- dist.taxo(data[[s]])
  # Now compute Moran's I:
  n <- length(d)
  l <- p <- i <- vector(length = n)
  for(j in 1:n) {
    if(j == 1) Mat <- d[[j]]
    else Mat <- d[[j]] & !d[[j-1]]
    I.M  <- Moran.I(y, Mat, scale=TRUE);
    i[j] <- I.M$obs
    p[j] <- I.M$p.v
    l[j] <- names(d)[j]
  }

  # Create an object of class 'correlogram':
  corr <- list(obs=i, p.values=p, labels=l)
  class(corr) <- "correlogram"
  return(corr)
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
  if (class(phy) != "phylo") stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) stop("tree \" phy\" must have branch lengths.") 
  #Get the minimum and maximum distance in the tree:
  dist <- dist.phylo(phy)
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
  if (class(x) != "correlogram") stop("object \"x\" is not of class \"correlogram\"")
   # Draw the correlogram (using lattice library):
  library(lattice)
  # New Black & White device:
  # Black circles are significant at the 5% level:
  pch <- ifelse(x$p.values < test.level, 19, 21)
  # Plot it!
  return(xyplot(x$obs~ordered(x$l,levels=x$l), type="b", xlab="Rank", ylab="I / Imax", lty=2, lwd=2, cex=1.5, pch=pch, ...))
}

#co <- correlogram.formula(log10(SW) ~ Order/SuperFamily/Family/Genus, data=carn, test.level=0.01, ylim=c(-0.5,1))
#plot(co)


