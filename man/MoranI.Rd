\name{Moran.I}
\alias{Moran.I}
\title{Moran's I Autocorrelation Index}
\usage{
  Moran.I(x, dist, scaled = FALSE, na.rm = FALSE)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{dist}{a distance weigths matrix.}
  \item{scaled}{logical, indicate whether the index should be scaled to
    allow comparisons between indices (default to \code{FALSE}).}
	\item{na.rm}{a logical value indicating whether 'NA' values in 'x'
		should be stripped before the computation proceeds.}
}
\description{
  Computes Moran's I autocorrelation index of a given vector \code{x}
  according to a distance weights matrix \code{dist} using the method
  described by Gittleman and Kot (1990).
}
\details{
  The \code{dist} matrix is used as distance weights, and Moran's I indice
	is computed using the formula:
  \deqn{
		I = \frac{n}{S_0} \frac{\sum_{i=1}^n\sum_{j=1}^n w_{i,j}(y_i - \overline{y})(y_j - \overline{y})}{\sum_{i=1}^n {(y_i - \overline{y})}^2}
	}{
		\code{I = n/S0 * (sum\{i=1..n\} sum\{j=1..n\} wij(yi - ym))(yj - ym) / (sum\{i=1..n\} (yi - ym)^2)}
	}
  with
	\itemize{
		\item \eqn{y_i}{yi} = observations
		\item \eqn{w_{i,j}}{wij} = distance weight
		\item \eqn{n} = number of observations
		\item \eqn{S_0}{S0} = \eqn{\sum_{i=1}^n\sum_{j=1}^n wij}{\code{sum_{i=1..n} sum{j=1..n} wij}}
	}
  The value of I may be evaluated under two null hypotheses, normality and randomization.
  Only the randomization hypothesis is implemented for now.

	The 'na.rm' option is forwarded to all 'sd', 'sum' and 'mean' functions used in computations.
	Note that 'NA' values are allowed only in 'x' and not in 'dist'.
}
\value{
  A list containing the attributes:
  \item{observed}{Moran's I index of \code{x}.}
  \item{expected}{Expected value of I under the null hypothesis.}
  \item{sd}{The standard deviation of I under the null hypothesis.}
  \item{p.value}{The p-value of having the observed value under the null hypothesis.}
}
\references{
  Gittleman, J. L. and Kot, M. (1990) Adaptation: statistics and a null
  model for estimating phylogenetic effects. \emph{Systematic Zoology},
  \bold{39}, 227--241.
}
\author{Julien Dutheil \email{julien.dutheil@univ-montp2.fr}}
\seealso{
	\code{\link{cophenetic.phylo}}, \code{\link{discrete.dist}}, \code{\link{weight.taxo}}
}
\examples{
### This will draw the correlogram from Gittleman and Kot 1990
### This is only for example, consider function plot.correlogram()
### to draw correlograms.
data(carnivora)
attach(carnivora)
# Compute distance matrix for each taxonomic level:
dG  <- weight.taxo(Genus)
dF  <- weight.taxo(Family)
dSF <- weight.taxo(SuperFamily)
dO  <- weight.taxo(Order)
# We draw the correlogram of average body weights,
# with log10-transformed variable:
IG  <- Moran.I(log10(SW), dG        , scale=TRUE)
IF  <- Moran.I(log10(SW), dF  & !dG , scale=TRUE)
ISF <- Moran.I(log10(SW), dSF & !dF , scale=TRUE)
IO  <- Moran.I(log10(SW), dO  & !dSF, scale=TRUE)
# All Moran's I indices:
i <- c(IG$obs, IF$obs, ISF$obs, IO$obs)
# With their corresponding p-values:
p <- c(IG$p.v, IF$p.v, ISF$p.v, IO$p.v)
# Here's the legend:
l <- c("G", "F", "SF", "O")
# Draw the correlogram (using lattice library):
library(lattice)
# New Black & White device:
trellis.device(device=X11, new=FALSE, color=FALSE)
# Black circles are significant at the 5% level:
pch <- ifelse(p < 0.05, 19, 21)
# Plot it!
xyplot(i~ordered(l,levels=l), type="b", xlab="Rank", ylab="I / Imax",
    lty=2, lwd=2, cex=1.5, pch=pch, ylim=c(-0.75,0.75))
}
\keyword{models}
\keyword{regression}
