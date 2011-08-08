\name{mantel.test}
\alias{mantel.test}
\title{Mantel Test for Similarity of Two Matrices}
\usage{
mantel.test(m1, m2, nperm = 1000, graph = FALSE,
            alternative = "two.sided",  ...)
}
\arguments{
  \item{m1}{a numeric matrix giving a measure of pairwise distances,
    correlations, or similarities among observations.}
  \item{m2}{a second numeric matrix giving another measure of pairwise
    distances, correlations, or similarities among observations.}
  \item{nperm}{the number of times to permute the data.}
  \item{graph}{a logical indicating whether to produce a summary graph
    (by default the graph is not plotted).}
  \item{alternative}{a character string defining the alternative
    hypothesis:  \code{"two.sided"} (default),  \code{"less"},
    \code{"greater"}, or any unambiguous abbreviation of these.}
  \item{\dots}{further arguments to be passed to \code{plot()} (to add a
    title, change the axis labels, and so on).}
}
\description{
  This function computes Mantel's permutation test for similarity of two
  matrices. It permutes the rows and columns of the two matrices
  randomly and calculates a \eqn{Z}-statistic.
}
\details{
  The function calculates a \eqn{Z}-statistic for the Mantel test, equal to
  the sum of the  pairwise product of the lower triangles of the
  permuted matrices, for each permutation of rows and columns. It
  compares the permuted distribution with the \eqn{Z}-statistic observed
  for the actual data.

  If \code{graph = TRUE}, the functions plots the density estimate of
  the permutation distribution along with the observed \eqn{Z}-statistic
  as a vertical line.

  The \code{\dots} argument allows the user to give further options to
  the \code{plot} function: the title main be changed with \code{main=},
  the axis labels with \code{xlab =}, and \code{ylab =}, and so on.
}
\value{
  \item{z.stat}{the \eqn{Z}-statistic (sum of rows*columns of lower
    triangle) of the data matrices.}
  \item{p}{\eqn{P}-value (quantile of the observed \eqn{Z}-statistic in
    the permutation distribution).}
  \item{alternative}{the alternative hypothesis.}
}
\references{
  Mantel, N. (1967) The detection of disease clustering and a
  generalized regression approach. \emph{Cancer Research}, \bold{27},
  209--220.

  Manly, B. F. J. (1986) \emph{Multivariate statistical methods: a primer.}
  London: Chapman & Hall.
}
\author{
  Original code in S by Ben Bolker \email{bolker@zoo.ufl.edu}, ported
  to \R by Julien Claude \email{claude@isem.univ-montp2.fr}
}
\examples{
q1 <- matrix(runif(36), nrow = 6)
q2 <- matrix(runif(36), nrow = 6)
mantel.test(q1, q2, graph = TRUE,
            main = "Mantel test: a random example with 6 X 6 matrices",
            xlab = "z-statistic", ylab = "Density",
            sub = "The vertical line shows the observed z-statistic")
}
\keyword{multivariate}
