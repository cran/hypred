\name{hypredSNPdist}
\docType{methods}
\alias{hypredSNPdist}
\alias{hypredSNPdist,hypredGenome-method}

\title{
Find the genetic distance and recombination frequency between two SNP on
a chromosome
}

\description{
A generic function to find the genetic distance and recombination frequency between two SNP on
a chromosome.
}

\usage{
hypredSNPdist(object, \dots)
\S4method{hypredSNPdist}{hypredGenome}(object, chromosome, SNP1, SNP2,
block = FALSE)
}
\arguments{
  \item{object}{an object of a class that holds information on map
    positions of loci, typically an \code{"\linkS4class{hypredGenome}"} object.}
  \item{chromosome}{The chromosome on which the SNP are, given as integer.}
  \item{SNP1}{The first SNP, given as integer.}
  \item{SNP2}{The second SNP, given as integer. This SNP must have a
    higher index than \code{SNP1}.}
  \item{block}{logical argument indicating if the map positions of the
    SNP themselfes should be used (\code{FALSE}) or the map positions of the block in
    which they are located (\code{TRUE}).}
  \item{\dots}{Methods may require further arguments.}
}
\value{
  A two row, one column matrix with the genetic distance (in M) in row
  one and the recombination frequency in row two.
}
\details{
The recombination frequency is calculated by the inverse of the Haldane
mapping function.
}

\seealso{
The function \code{\link{hypredDInv}} which is used by the
\code{hypredGenome} method.
}


\examples{

## one chromosome of length 1 M and 5 SNP

map <- hypredGenome(1, 1.0, 5)

## distance and recombination frequency between SNP 1 and 3 on
## chromosome 1
hypredSNPdist(map,1,1,3)

}

\keyword{methods}
\keyword{univar}


