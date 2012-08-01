\name{hypredNewMap}
\docType{methods}
\alias{hypredNewMap}
\alias{hypredNewMap,hypredGenome-method}

\title{
Modify the map of a hypredGenome object
}

\description{
A generic function to modify the genetic map defined in the
\code{"\linkS4class{hypredGenome}"} object
}

\usage{
hypredNewMap(object, \dots)
\S4method{hypredNewMap}{hypredGenome}(object, new.map)
}
\arguments{
  \item{object}{an object of a class that holds information on map
    positions of loci, typically an \code{"\linkS4class{hypredGenome}"} object.}
  \item{new.map}{numeric vector with the map positions of all loci
    (marker and qtl)}
  \item{\dots}{Methods may require further arguments.}
}
\value{
A \code{"\linkS4class{hypredGenome}"} object, modified accordingly.
}
\details{
Only the map positions can be changed, the length of the chromosomes,
their number and the number of SNP and QTL per chromosome is not
modified. There are some validity checks included to assure that the new
map agrees with these parameters.

The vector given to \code{new.map} holds the positions of all loci. So
if, for example, the number of chromosomes is three and the number of
loci per chromosome is 10, then \code{new.map[1:10]} must contain the
map of chromosome 1, \code{new.map[11:20]} the map of chromosome 2 and
\code{new.map[21:30]} the mao of chromosome 3.

The postions must give the distance to the beginning of the chromosome,
not the distance between adjacent loci!

If the object included SNP blocks, they will not be adapted! Instead all
slots referring to the blocks are set to \code{NULL} and a warning is issued. 
}

\seealso{
The function \code{\link{hypredGenome}} which creates the original object.
}


\examples{

## one chromosome of length 1 M and 5 SNP

random_map <- hypredGenome(1, 1.0, 5)

## map positions are random

slot(random_map,"pos.snp")

## change to regular positions

change_map <- c(0.1, 0.2, 0.3, 0.4, 0.5)

regular_map <- hypredNewMap(random_map, change_map)

slot(regular_map,"pos.snp")

}

\keyword{methods}

\keyword{datagen}


