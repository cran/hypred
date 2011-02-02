\name{hypredRecombine}
\docType{methods}
\alias{hypredRecombine}
\alias{meiosisFUNallChr}
\alias{hypredRecombine,hypredGenome-method}

\title{
  Recombine two haploid genomes
}

\description{
A generic function to recombine two haploid parental genomes (as found
in sex cell) into one haploid progeny genome (i.e. the gamete)
}

\usage{
hypredRecombine(object, \dots)
\S4method{hypredRecombine}{hypredGenome}(object, genomeA, genomeB,
mutate, mutation.rate.snp, mutation.rate.qtl,  block)
}
\arguments{
  \item{object}{an object of a class that holds information on genome
    parameters necessary to simulate meiosis, typically an
    \code{"\linkS4class{hypredGenome}"} object.}
  \item{genomeA}{One of the haploid parental genomes. An integer vector
    or one row matrix. The two SNP alleles are coded as 1 and 0.}
  \item{genomeB}{The other of the haploid parental genomes. An integer vector
    or one row matrix. The two SNP alleles are coded as 1 and 0.}
  \item{mutate}{Logical argument indicating if mutations should be simulated as well,
    \code{TRUE} or \code{FALSE}.}
  \item{mutation.rate.snp}{The mutation rate of SNP \eqn{\mu} (per base pair), a numeric
    value.  If \code{mutate = FALSE} this argument is ignored.}
  \item{mutation.rate.qtl}{The mutation rate of QTL \eqn{\mu} (per base pair), a numeric
    value.  If \code{mutate = FALSE} this argument is ignored.}
  \item{block}{logical argument indicating if the map positions of the
    SNP themselfes should be used (\code{FALSE}) in the simulation or
     the map positions of the block in which they are located
     (\code{TRUE}).}
   \item{\dots}{Methods may require further arguments.}
 }
 
 \value{
   A one row matrix with \code{ncol(genomeA)} resp. \code{ncol(genomeB)}
   columns, giving the genotype of the gamete.  The two SNP alleles are coded as 1 and 0.
 }

 \details{ The meiosis is simulated according to the count-location
   model that assumes no interference. Here the number of crossing over on 
   a random meiotic product of a chromosome is Poission distributed with
   parameter \eqn{\lambda = L} (\eqn{L} being the length of the
   chromosome in Morgan). The location of the crossing over follow a
   uniform distribution over the interval \eqn{(0,L)}.

   The computations for simulating the meiosis are performed by a
   routine written in C and called through \code{.C}.

   Since the count-location model corresponds to the Haldane mapping
   function, observed recombination frequencies between loci must agree
   with the expected ones as calculated from the inverse of the Haldane
   function from the map distances. This fact was used to validate the
   algorithm. If you want to do this yourself, remember to use a
   sufficently high number of progenies.
   
   The occurence of mutations is simulated according to a simple
   count-location model as well. The number of observed SNP mutations in
   the gamete follows a Poission distribution with parameter
   \eqn{\lambda =} \eqn{\mu} times the number of SNP (i.e. loci not
   assigned QTL) on a chromosome times the number of chromosomes
   (\code{mutation.rate.snp * object@num.snp.chr *
   object@num.chr}). Then each SNP-loci has an equal chance of beeing
   mutated, irrespective on whether mutations occured in its vicinity. A
   mutation event always changes a 1 allele to a 0 allele and vice versa
   (i.e. 0 to 1; 1 to 0). During one meiosis, a given SNP-loci can only
   be mutated twice if the number of mutations exceeds the total number
   of SNP. This should rarely happen and is included only to avoid an
   error. The same applies to mutations of QTL.

   Since all the simulated loci are treated as SNP, \eqn{\mu} should be
   interpreted as \eqn{\mu} per \emph{base-pair} per meiosis.

   

 }

\references{
Wu, R.; Chang-Xing, M; Casella, G.; (2007) \emph{Statistical Genetics of
Quantitative Traits}.
New York: Springer.
}

\author{
Frank Technow
}

\note{ The function call the C routine \code{meiosisFUNallChr} to
  simulate the meiosis.}

\seealso{
  The function \code{\link{hypredGenome}} which creates the
  \code{"\linkS4class{hypredGenome}"} object, and the function
  \code{\link{hypredFounder}} which creates genomes of founder lines
  from which progenies can be created by function \code{hypredRecombine}.
}


\examples{

## one chromosome of length 1 M and 5 SNP

genomeDef <- hypredGenome(1, 1.0, 5)

## produce two haploid founder line genomes

founder <- hypredFounder(genomeDef,1)
founder

## produce one gamete from a cross between the two
set.seed(135)
gamete <- hypredRecombine(genomeDef,
                          genomeA = founder[1,],
                          genomeB = founder[2,],
                          mutate = FALSE,
                          block = FALSE)

print(gamete)
}

\keyword{methods}
\keyword{datagen}
