\name{hypredCode}
\docType{methods}
\alias{hypredCode}
\alias{hypredCode_FUN_dom}
\alias{hypredCode,hypredGenome-method}

\title{
  Create various genomic design matrices
}

\description{ A generic function to create various kinds of genomic
design matrices, for estimation of substitution effects, additive and
dominance effects or genome specific effects of SNPs}


\usage{
hypredCode(object, \dots)
\S4method{hypredCode}{hypredGenome}(object, genotypes, DH, type)
}
\arguments{
  \item{object}{an object of a class that holds information on genome
    parameters necessary to create design matrices, typically an
    \code{"\linkS4class{hypredGenome}"} object.}
  
  \item{genotypes}{ Integer matrix giving the genotypes of the
    individuals (as a series of 1s and 0s). If \code{DH = FALSE}, always
    two adjacent rows stand for one individual (the two chromosome
    sets).  If \code{DH = TRUE}, always one row stands for one
    individual, since then the chromosome sets are identical copies of
    each other.}

  \item{DH}{logical argument indicating if one row in the matrix given to \code{genotypes}
    belongs to one  DH line, i.e. completely homozygous line
    (\code{TRUE}) or if two adjacent rows belong to one, DH or not, individual (\code{FALSE}).}

  \item{type}{character string giving the type of design matrix to
    create \code{c("012", "-101", "Xu2003", "genome.specific")}, see the
    \code{Details} section. \code{"genome.specific"} can only be used
    when all odd rows in the matrix given to \code{genotypes} belong to
    genome 1 and all even rows to genome 2!}
  \item{\dots}{Methods may require further arguments.}
}
 
 \value{The desired design matrix, whose dimensions depend on
 \code{type}. The column names of the matrix correspond to the loci IDs}

 \details{
   The SNP loci that represent QTL are not included in the design
   matrices, exept when they are perfect markers.

   
   Description of different design matrices, (\eqn{N} = number of
   individuals (observations), \eqn{M} = number of markers),
   (sub)matrices that code additive effects are denoted as \eqn{X},
   (sub)matrices that code dominance effects are denoted as \eqn{W}. 
   \eqn{i} indexes the individual, \eqn{j} the marker locus.

   \describe{
     
     \item{\code{type = "012"}}{Standard design matrix for estimation of
       substitution effects. \eqn{X} has dimensions [\eqn{N}, \eqn{M}]. \eqn{x_{ij} = 2}
       if the locus is 11, \eqn{x_{ij} = 1} if the locus is 10(01) and
       \eqn{x_{ij} = 0} if the locus is 00.}
       
     \item{\code{type = "-101"}}{Design matrix for estimation of
       substitution effects. \eqn{X} has dimensions [\eqn{N}, \eqn{M}].\eqn{x_{ij} = 1} if
       the locus is 11, \eqn{x_{ij} = 0} if the locus is 10(01) and
       \eqn{x_{ij} = -1} if the locus is 00.}
     
     \item{\code{type = "Xu2003"}}{Design matrix for estimation of
       additive and dominance effects. The matrix is partitioned as
       \eqn{[\eqn{X} \eqn{W}]}. The dimensions of \eqn{X} and \eqn{W} are [\eqn{N}, \eqn{M}].
       \eqn{x_{ij} = \sqrt{2}} and \eqn{w_{ij} = -1} if the locus is 11,
       \eqn{x_{ij} = 0} and \eqn{w_{ij} = 1} if the locus is 10(01),
       \eqn{x_{ij} = -\sqrt{2}} and \eqn{w_{ij} = -1} if the locus is
       00. This is the coding scheme used by Xu (2003).}

     \item{\code{type = "genome.specific"}}{Design matrix for estimation
       of genome specific additive and dominance effects. The matrix is
       partitioned as \eqn{[\eqn{X(1)} \eqn{X(2)} \eqn{W(1)} \eqn{W}(2)]}. The dimensions of
       each submatrix are [\eqn{N}, \eqn{M}]. \eqn{X(1)} and \eqn{W(1)} code the
       additive and dominance effects of genome 1, \eqn{X(2)} and
       \eqn{W(2)} code the additive and dominance effects of genome
       2. An element of \eqn{X(1)} or \eqn{X(2)} is 1 if a 1 allele is
       contributed by the respective genome, else it is
       zero. \eqn{x(1)_{ij} = 1} if a 1 allele at the locus is
       contributed by genome 1, \eqn{x(1)_{ij} = 0} if not, same for
       genome 2. \eqn{w(1)_{ij} = 1} if the locus is heterozygous and
       the 1 allele comes from genome 1, else \eqn{w(1)_{ij} = 0}, same
       for genome 2.

       Example:
       \describe{
	 \item{genotype 10}{ \eqn{x(1)_{ij} = 1},
	   \eqn{x(2)_{ij} = 0}, \eqn{w(1)_{ij} = 1} and \eqn{w(2)_{ij} = 0}}

	 \item{genotype 01}{ for the genotype 01, \eqn{x(1)_{ij} = 0},
	   \eqn{x(2)_{ij} = 1}, \eqn{w(1)_{ij} = 0} and \eqn{w(2)_{ij} = 1}}

	 \item{genotype 11}{ for the genotype 11, \eqn{x(1)_{ij} = 1},
	   \eqn{x(2)_{ij} = 1}, \eqn{w(1)_{ij} = 0} and \eqn{w(2)_{ij} = 0}}

	 \item{genotype 00}{ for the genotype 00, \eqn{x(1)_{ij} = 0},
	   \eqn{x(2)_{ij} = 0}, \eqn{w(1)_{ij} = 0} and \eqn{w(2)_{ij} = 0}}
       }

       \emph{For this to work, all odd rows in the matrix given to
       \code{genotypes} must belong to genome 1 and all even rows to
       genome 2. If this is not the case, results will be nonsense!}
     
    }
  }
}

\author{
  Frank Technow
}

\note{The function call the C routine \code{hypredCode_FUN_dom} to code
  the "Xu2003" type matrix.}

\seealso{ The function \code{\link{hypredRecombine}} which is used to
  create progeny genomes, the function \code{\link{hypredNewQTL}} which
  allows to assign new QTLs to the \code{"\linkS4class{hypredGenome}"}
  object and the function \code{\link{hypredTruePerformance}} which
  determines the standart, non-genome-specific values.}


\examples{

## one chromosome of length 1 M and 5 SNP
genomeDef <- hypredGenome(1, 1.0, 5)

## assign one QTL with and additive effect of 1 to the first locus

genomeDef <- hypredNewQTL(genomeDef,
                          new.id.add = 1,
                          new.eff.add = 1)

summary(genomeDef)


## produce two haploid founder line genomes
founder <- hypredFounder(genomeDef,1)
founder

## produce two progeny from a cross between the two
## (this corresponds to F2 individuals)

set.seed(134)

gamete1 <- hypredRecombine(genomeDef,
                          genomeA = founder[1,],
                          genomeB = founder[2,],
                          mutate = FALSE,
                          block = FALSE)

gamete2 <- hypredRecombine(genomeDef,
                          genomeA = founder[1,],
                          genomeB = founder[2,],
                          mutate = FALSE,
                          block = FALSE)

gamete3 <- hypredRecombine(genomeDef,
                          genomeA = founder[1,],
                          genomeB = founder[2,],
                          mutate = FALSE,
                          block = FALSE)

gamete4 <- hypredRecombine(genomeDef,
                          genomeA = founder[1,],
                          genomeB = founder[2,],
                          mutate = FALSE,
                          block = FALSE)


F2 <- rbind(gamete1, gamete2, gamete3, gamete4)

print(F2)

## type "012" matrix

mat012 <- hypredCode(genomeDef,
                     genotypes = F2,
                     DH = FALSE,
                     type = "012")

## note that loci 1 is not included
print(mat012)


## type "genome.specific" matrix

matGS <- hypredCode(genomeDef,
                     genotypes = F2,
                     DH = FALSE,
                     type = "genome.specific")

## note that loci 1 is not included
print(matGS)

}

\keyword{methods}
\keyword{datagen}
