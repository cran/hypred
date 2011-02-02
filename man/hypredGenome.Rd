\name{hypredGenome}
\alias{hypredGenome}
\title{Define genome parameters}

\description{ This function is used to define the genome parameters:
number and length of chromosomes, number of SNP markers and haplotype
blocks. This function also performes part of the simulation, by
simulating a genetic map (which can be modified later).}


\usage{ hypredGenome(num.chr, len.chr, num.snp.chr, num.blocks.chr = NULL) }

\arguments{
  \item{num.chr}{
    Integer argument indicating the number of chromosomes.
  }
  \item{len.chr}{
    Numeric vector giving the length of each chromosome in Morgan. This
    vector must have as many elements as there are chromosomes. 
  }
  \item{num.snp.chr}{
    Integer value indicating the number of SNP markers per chromosome. 
  }

  \item{num.blocks.chr}{
    Integer value indicating the number of haplotype blocks per
    chromosome. Default is \code{NULL}.
  }
}

\details{
  \describe{   
    \item{loci}{The markers are interpreted and treated as SNP with two
      alleles, 1 and 0.  The number of SNP marker is the same for each
      chromosome.}
    
    \item{map positions}{The map positions of all loci are sampled from
      the uniform distribution over the interval \eqn{(0,L)}, where \eqn{L}
      is the length of the chromosome in Morgan.}
    \item{haplotype blocks}{
	All SNP that are located on one block are treated as
	perfectly linked. All functions that use the map positions have
	arguments to indicate wheter the block positions of QTL and SNP
	should be used or their individual positions.

	If the argument is non \code{NULL}, then the following condition must hold, \cr
	\code{(num.blocks.chr/num.snp.chr) < 0.1001}, else
	there is an error.

	The number of loci on all blocks but the last is:\cr
	\code{floor(num.snp.chr/num.block.chr)}\cr
	If \code{num.snp.chr} mod \code{num.block.chr} != 0, 
	the last block will contain more loci than the other blocks.
	
	The map position of a block is the one of its
	center SNP. These center SNP are found as follows:
    
	block1: \code{floor(num.snp.chr/
	  (2*num.block.chr))}
    
	block2: \code{floor(num.snp.chr/
	  (2*num.block.chr))} + \cr
	\code{1* floor(num.snp.chr/num.block.chr)}

	block3: \code{floor(num.snp.chr/
	  (2*num.block.chr))} + \cr
	\code{2 * floor(num.snp.chr/num.block.chr)}
	
	and so on.
      }
      }
}

\value{
  An object of class \code{"\linkS4class{hypredGenome}"}
}

\author{
Frank Technow
}

\note{
  Some of the results of this function are subject to randomness. To obtain
  reproducible results, you should use \code{set.seed()} to set a seed
  for the random number generator.
}


\seealso{
  The function \code{\link{hypredNewMap}} which allows to modify the
  genetic map. The function \code{\link{hypredNewQTL}} which allows to assign
  QTL or modify already existing ones.
}


\examples{
gen.def <- hypredGenome(num.chr = 2, 
                        len.chr = c(1.0,1.2),
                        num.snp.chr = 90)


summary(gen.def)

}
\keyword{datagen}

