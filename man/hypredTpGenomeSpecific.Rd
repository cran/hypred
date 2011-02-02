\name{hypredTpGenomeSpecific}
\docType{methods}
\alias{hypredTpGenomeSpecific}
\alias{hypredTpGenomeSpecific,hypredGenome-method}

\title{
 Determine the true genetic performance (genotypic value) of individuals 
}

\description{ A generic function to determine the true genetic
performance (genotypic value) of individuals according to their genome,
with effects beeing genome specific.}


\usage{
hypredTpGenomeSpecific(object, \dots)
\S4method{hypredTpGenomeSpecific}{hypredGenome}(object, genotypes,
add.eff.G1 = NULL, add.eff.G2 = NULL, dom.eff.G1 = NULL, dom.eff.G2 = NULL)
}
\arguments{
  \item{object}{an object of a class that holds information on genome
    parameters necessary to determine the true genetic performance of individuals, typically an
    \code{"\linkS4class{hypredGenome}"} object.}
  
  \item{genotypes}{Integer matrix giving the genotypes of the
    individuals (as a series of 1s and 0s). Two adjacent rows stand for
    one individual (the two chromosome sets). Consistently, odd rows
    belong to genome 1, even rows to genome 2}
  
  \item{add.eff.G1}{numeric vector with the additive effects of the QTL
    alleles from genome 1.}

  \item{add.eff.G2}{numeric vector with the additive effects of the QTL
    alleles from genome 2.}

  \item{dom.eff.G1}{numeric vector with the dominance effects of the
    QTL, when the loci is heterozygous and the 1 allele comes from genome 1.}

  \item{dom.eff.G2}{numeric vector with the dominance effects of the
    QTL, when the loci is heterozygous and the 1 allele comes from
    genome 2.}
  
  \item{\dots}{Methods may require further arguments.}
 }
 
 \value{ A one column matrix with as many elements as there were
   individuals, giving their true genetic performances in the order in
   which they appeard in the the matrix.}

 \details{ The additive contribution of a genome to the genotypic value
   is determined by summing the additive effects of this genome of all
   loci where the 1 allele comes from that genome.

   The contribution of a genome to the dominance part of the genotypic
   value is determined by summing the dominance effects of all dominance
   QTL that are heterozygous \emph{and} where the 1 allele comes from
   that genome.

   Example: Assume that there is only one QTL, with both additive and
   dominance effect. The additive effect is 1 for genome 1 and 0.5 for
   genome 2, the dominance effect is 0.5 for genome 1 and -0.5 for
   genome 2. An individual with genotype 0-1 at this QTL would have a
   performance of 0 (additive part genome 1: 0, additive part genome 2:
   0.5, dominance part genome 1: 0, dominance part genome 2: -0.5 ), an
   individual with 1-0 would have a performance of 1.5 (additive part
   genome 1: 1, additive part genome 2: 0, dominance part genome 1: 0.5,
   dominance part genome 2: 0).

   The defaults of all four arguments that give the effects are
   \code{NULL}. When an argument is \code{NULL}, the effects are taken
   from the object given to \code{object}. Hence, calling
   \code{hypredTpGenomeSpecific} with all effect arguments \code{NULL},
   produces identical results as \code{hypredTruePerformance}.

   An indirect way of obtaining the additive genotypicic values when
   dominance QTL were assigned is to call the function with
   \code{add.eff.G1 = NULL}, \code{add.eff.G2 = NULL}, \code{dom.eff.G1
   = rep(0, number.dominance.qtl)}, \code{dom.eff.G2 = rep(0,
   number.dominance.qtl)}

 The number of effects must agree with the numbers of QTL assigned to
 \code{object}. If no QTLs are assigned, the function returns an error.  }

\author{
Frank Technow
}
 
\seealso{ The function \code{\link{hypredRecombine}} which is used to
  create progeny genomes, the function \code{\link{hypredNewQTL}} which
  allows to assign new QTLs to the \code{"\linkS4class{hypredGenome}"}
  object and the function \code{\link{hypredTruePerformance}} which
  determines the standart, non-genome-specific values.}


\examples{

## one chromosome of length 1 M and 5 SNP

genomeDef <- hypredGenome(1, 1.0, 5)

genomeDef <- hypredNewQTL(genomeDef,
                          new.id.add = 1,
                          new.eff.add = 1,
                          new.id.dom = 1,
                          new.eff.dom = 0.5)


summary(genomeDef)

## assign one QTL with and additive effect of 1 to the first loci,
## with the also having a dominance effect of 0.5 . Firts QTL
## is perfect marker as well.

genomeDef <- hypredNewQTL(genomeDef,
                          new.id.add = 1,
                          new.id.dom = 1,
                          new.id.per.mar = 1,
                          new.eff.add = 1,
                          new.eff.dom = 0.5)

summary(genomeDef)

## QTL genotype individual 1 = 1-1
## QTL genotype individual 2 = 0-1

individuals <- matrix(c(1,0,1,0,1,
                        1,1,0,1,0,
                        0,0,1,0,1,
                        1,1,0,1,0), nrow = 4, byrow = TRUE)

gvalues <- hypredTpGenomeSpecific(genomeDef,
                                  genotypes = individuals,
                                  add.eff.G1 = 1, ## additive effect from genome 1
                                  add.eff.G2 = 0.5, ## additive effect from genome 2
                                  dom.eff.G1 = 0.5, ## dominance effect from genome 1 
                                  dom.eff.G2 = -0.5) ## dominance effect from genome 2 

## value of ind. 1 must be 1.5, of ind. 2 must be 0
print(gvalues)

## test
gvalues <- as.vector(gvalues)

stopifnot(all.equal(gvalues, c(1.5, 0.0)))

}

\keyword{methods}
\keyword{datagen}
