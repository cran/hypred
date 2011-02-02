\name{hypredFounder}
\docType{methods}
\alias{hypredFounder}
\alias{hypredFounder,hypredGenome-method}

\title{
  Creates gametes of two founder lines
}

\description{ A generic function to create the gametes of two founder
lines that differ at each locus, according to a defined genome.}

\usage{
hypredFounder(object, \dots)
\S4method{hypredFounder}{hypredGenome}(object, prob.snp)
}
\arguments{
  \item{object}{an object of a class that holds information on a genome
    typically an \code{"\linkS4class{hypredGenome}"} object.}
  
  \item{prob.snp}{Numeric value that gives the probability that founder
    line 1 has the 1 allele at a locus.}
  \item{\dots}{Methods may require further arguments.}
}

\value{
  A two row matrix with as many columns as there are loci. Row 1 holds
  the gamete of line 1, row 2 the one of line 2. Row 2 is the exact
  negative of row 1 (i.e. if \code{value[1, 1] == 1}, than \code{value[2, 1] == 0}).
  }
  
  \details{ For each locus a 1 is sampled (using function \code{sample})
    with probability \code{prob.snp}.  Setting \code{prob.snp = 1},
    results in gamete 1 having all the one alleles.
  }

  \seealso{ The function \code{\link{hypredGenome}} which creates the
    object, the function \code{\link{hypredNewMap}} which allows to
    modify the genetic map and the function \code{\link{hypredNewQTL}}
    which allows to assign QTL or modify already existing ones. The
    genome of the founders can be a starting point for the simulation of
    recombination using function \code{\link{hypredRecombine}}.}
  


\examples{

## one chromosome of length 1 M and 5 SNP

genome <- hypredGenome(1, 1.0, 5)

founders <- hypredFounder(genome,
                          prob.snp = 0.5) 

print(founders)

}

\keyword{methods}
\keyword{datagen}


