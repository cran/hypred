\name{hypredTruePerformance}
\docType{methods}
\alias{hypredTruePerformance}
\alias{hypredTruePerformance,hypredGenome-method}

\title{
 Determine the true genetic performance (genotypic value) of individuals 
}

\description{ A generic function to determine the true genetic
performance (genotypic value) of individuals according to their genome.}


\usage{
hypredTruePerformance(object, \dots)
\S4method{hypredTruePerformance}{hypredGenome}(object, genotypes, DH)
}
\arguments{
  \item{object}{an object of a class that holds information on genome
    parameters necessary to determine the true genetic performance of individuals, typically an
    \code{"\linkS4class{hypredGenome}"} object.}
  \item{genotypes}{Integer matrix giving the genotypes of the
    individuals (as a series of 1s and 0s). If \code{DH = FALSE}, always
    two adjacent rows stand for one individual (the two chromosome sets).
    If \code{DH = TRUE}, always one row stands for one individual, since
    then the chromosome sets are identical copies of each other. }
  \item{DH}{logical argument indicating if one row in the matrix given to \code{genotypes}
    belongs to one  DH line, i.e. completely homozygous line
    (\code{TRUE}) or if two adjacent rows belong to one, DH or not,
    individual (\code{FALSE}).}
  \item{\dots}{Methods may require further arguments.}
 }
 
 \value{
   A one column matrix with as many elements as there were individuals,
   giving their true genetic performances in the order in which they
   appeard in the the matrix.  
 }

 \details{ The additive genotypic value is determined by multiplying the
   effects of the QTL with the number of 1 alleles of the given QTL in
   an individual and then summing over the QTL. The part of the
   genotypic value due to dominance is determined by summing the
   dominance effects of all dominance QTL that are heterozygous.

   Example: Assume that there is only one QTL, with both additive and
   dominance effect. The additive effect is 1, the dominance effect is
   0.5. An individual with genotype 0-1 at this QTL would have a
   performance of 1.5, a individual with 0-0 a performance of 0 and an
   individual with 1-1 a performance of 2.

   Making use of \code{DH = TRUE}, when the individuals under
   consideration are indeed DHs, has two advantages, one is that the
   comupations run a littly faster, because dominance can be
   ignored. The main advantage, however, is that only half the amount of
   data needs to be carried around because one row in the martrix is
   sufficient to represent the genotype. When two adjacent (identical)
   rows represent the genotypes of the DH individuals, using \code{DH =
   FALSE} gives identical results to \code{DH = TRUE} when one row
   stands for one DH individual.

   \emph{When the individuals are not DH (or not all of them), using \code{DH
   = TRUE} would create nonsense results!}
   
 }

\author{
Frank Technow
}
 
\seealso{ The function \code{\link{hypredRecombine}} which is used to
  create progeny genomes and the function \code{\link{hypredNewQTL}}
  which allows to assign new QTLs to the
  \code{"\linkS4class{hypredGenome}"} object.}


\examples{

## one chromosome of length 1 M and 5 SNP
genomeDef <- hypredGenome(1, 1.0, 5)

## assign one QTL with and additive effect of 1 to the first two loci,
## with the first loci also having a dominance effect of 0.5 . Firts QTL
## is perfect marker.

genomeDef <- hypredNewQTL(genomeDef,
                          new.id.add = 1:2,
                          new.id.dom = 1,
                          new.id.per.mar = 1,
                          new.eff.add = c(1,1),
                          new.eff.dom = 0.5)

summary(genomeDef)

## QTL genotype individual 1 = 1-1; 0-1 value must be 3
## QTL genotype individual 2 = 1-0; 0-1 value must be 2.5

individuals <- matrix(c(1,0,1,0,1,
                        1,1,0,1,0,
                        1,0,1,0,1,
                        0,1,0,1,0), nrow = 4, byrow = TRUE)

gvalues <- hypredTruePerformance(genomeDef,
                                genotypes = individuals,
                                DH = FALSE)

print(gvalues)

## test
gvalues <- as.vector(gvalues)

stopifnot(all.equal(gvalues, c(3,2.5)))


}

\keyword{methods}
\keyword{datagen}
