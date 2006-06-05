\name{nj}
\alias{nj}
\title{Neighbor-Joining Tree Estimation}
\description{
  This function performs the neighbor-joining tree estimation of Saitou
  and Nei (1987).
}
\usage{
nj(X)
}
\arguments{
  \item{X}{a distance matrix; may be an object of class ``dist''.}
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new
  method for reconstructing phylogenetic trees. \emph{Molecular Biology
    and Evolution}, \bold{4}, 406--425.
}
\author{Emmanuel Paradis \email{paradis@isem.univ-montp2.fr}}
\seealso{
  \code{\link{write.tree}}, \code{\link{read.tree}},
  \code{\link{dist.dna}}
}
\examples{
### From Saitou and Nei (1987, Table 1):
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(NA, 8, 8)
M[row(M) > col(M)] <- x
M[row(M) < col(M)] <- x
diag(M) <- 0
rownames(M) <- colnames(M) <- 1:8
tr <- nj(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- nj(dist.dna(woodmouse))
plot(trw)
}
\keyword{models}
