\name{codeW.genomeSpecific}
\alias{codeW.genomeSpecific}
\title{
coding of genome specific design matrices for dominance effects
}

\description{ This function codes design matrix [W1,W2]. An element of
W1 will be 1 if the individual is heterozygous for this loci AND the SNP
comes from genome 1, else it is 0. An element of W2 will be 1 if the
individual is heterozygous for this loci AND the SNP comes from genome
2, else it is 0.}

\usage{
codeW.genomeSpecific(x, N, num.snp)
}
\arguments{ \item{x}{ integer matrix giving the genotypes of the
  individuals. Always two adjacent rows give the genotype of one
  individual, where the odd rows give genome 1 (from one of the parents)
  and even rows genome 2 (from the other parent). The presence of a SNP
  is coded as \code{as.integer{1}}, the absence as \code{as.integer{0}}.
  } 
  \item{N}{
    The number of individuals, integer scalar. This must be equal to \code{nrow(x)/2}
}
  \item{num.snp}{
    The number of SNP, integer scalar. This must be equal to \code{ncol(x)}
}
}

\value{
An integer matrix ([W1,W2]) with \code{N} rows and \code{num.snp * 2} columns. 
}

\author{
Frank Technow
}
\note{
This is an untility function and not intended to be
called directly by the user. 
}


\seealso{
 \code{\link{hypredCode}} and  \code{\link{hypredTpGenomeSpecific}}
 which use this function
}

\keyword{datagen}

