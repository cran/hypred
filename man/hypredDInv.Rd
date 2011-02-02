\name{hypredDInv}
\alias{hypredDInv}

\title{
Implements the inverse of the Haldane and Kosambi mapping functions
}
\description{
This function calculates recombination frequencies from genetic
distances according to the inverse of the Haldane or Kosambi mapping functions.
}
\usage{
hypredDInv(d, FUN)
}
\arguments{
  \item{d}{
    Genetic distance between loci in Morgan (M), \code{numeric} vector
}
  \item{FUN}{ Inverse function to be used -- \code{character} \code{"Haldane"}
    or \code{"Kosambi"}
}
}

\value{
numeric vector of \code{length(d)} giving the recombination frequencies
}
\references{
Wu, R.; Chang-Xing, M; Casella, G.; (2007) \emph{Statistical Genetics of
Quantitative Traits}.
New York: Springer.
}
\author{
Frank Technow
}


\seealso{
 \code{\link{hypredSNPdist}}
}
\examples{
## distance in Morgan between loci1 and loci2 and
## between loci2 and loci3

distances <- c(0.1,0.2)

hypredDInv(distances, "Haldane")
}
\keyword{univar}

