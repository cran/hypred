\name{hypredGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{hypredGenome-class}
\title{Class "hypredGenome"}
\description{
A class to hold information and data describing a genome in the context
of high density SNP data.
}
\section{Objects from the Class}{
  Objects are created by the function \code{hypredGenome}
  or by calls of the form \code{new("hypredGenome", ...)}
}
\section{Slots}{
  \describe{
    \item{\code{num.chr}:}{Object of class \code{"integer"}, giving the
      number of chromosomes. }
    \item{\code{len.chr}:}{Object of class \code{"numeric"}, giving the
      length of the chromosomes (a vector with one entry for each chromosome.) }
    \item{\code{num.snp.chr}:}{Object of class \code{"integer"}, giving
      the number of SNP marker per chromosome }
    \item{\code{pos.snp}:}{Object of class \code{"numeric"}, giving the
      map positions (in M) of each loci (SNP marker and QTL). This is
      one vector. The loci from chromosome 1 come first, then the ones
      from chromosome 2 and so on. If there are 100 loci per chromosome,
      then elements 1:100 give the map of chromosome 1, 101:200 the map of
      chromosome 2 and so on.}
    \item{\code{pos.snp.block}:}{Object of class \code{"numeric"} giving the
      map positions (in M) of each loci (SNP marker and QTL). This is
      one vector. The loci from chromosome 1 come first, then the one
      from chromosome 2 and so on. If there are 100 loci per chromosome,
      then elements 1:100 give the map of chromosome 1, 101:200 the map of
      chromosome 2 and so on. In contrast to \code{pos.snp}, not all
      loci have distinct map positions. The map positions of all loci
      within one block of loci are identical.}
    
    \item{\code{block.info}:}{Object of class \code{"list"} holding
      information about the mentioned blocks. There should be three
      elements in this list: \code{center}, an integer vector giving the
      loci that are in the center of the blocks; \code{end}, an
      integer vector giving the loci that are at the end of the
      blocks; \code{num.blocks.chr}, an integer value giving the number
      of such blocks on each chromosome.}

    \item{\code{chr.break.pts}:}{Object of class \code{"matrix"}. This
    is a two column matrix with as many rows as there are
    chromosomes that gives the ID of the first and last loci on a
    chromosome for vectors such as \code{pos.snp}. Chromosome 1 in row
    1, chromosome 2 in row 2 and so on.}
  
    \item{\code{num.add.qtl.chr}:}{Object of class \code{"integer"}. A
    single value that gives the number of QTL with additive effect on
    each chromosome (since all QTL must have an additve effect, this is
    identical to the number of QTL in general).}
  
  \item{\code{num.dom.qtl.chr}:}{Object of class \code{"integer"}. A
    single value that gives the number of QTL with a dominance effect on
    each chromosome (since all QTL must have an additve effect, this
    value must be equal or smaller than \code{num.add.qtl.chr}).}  
  
  \item{\code{pos.add.qtl}:}{Object of class \code{"list"}. This list
    holds information on the positions of the QTL with additive effects
    (hence all). There should be two elements in this list: \code{ID},
    an integer vector giving the IDs of the QTL; \code{M}, a
    numeric vector giving their map positions.}
    
  \item{\code{pos.dom.qtl}:}{Object of class \code{"list"}. Same as for
    \code{pos.add.qtl}, but only for QTL with dominance effects. Because
    only QTL with additive effects can have dominance effects, the
    elements must be subsets of the corresponding elements if \code{pos.add.qtl}.}

  \item{\code{add.and.dom.eff}:}{Object of class \code{"list"}.  There
    should be two elements in this list: \code{add}, a numeric vector with
    the additve effects; \code{dom}, a numeric vector with the dominance
    effects.}

  \item{\code{num.per.mar.chr}:}{Object of class \code{"integer"}. An
  integer value giving the number of perfect markers on each
  chromosome. A perfect marker is a QTL that is observed as a SNP.}

  \item{\code{per.mar.id}:}{Object of class \code{"integer"}. An integer
    vector giving the IDs of the perfect markers.}
}
}

\section{Methods}{
  \describe{
    \item{hypredCode}{\code{signature(object = "hypredGenome")}: Create
      design matrices of SNP effects for a set of individuals.}
    \item{hypredFounder}{\code{signature(object = "hypredGenome")}:
      Create the gametes of two homozygous founder lines that differ at each
      locus.}
    \item{hypredNewMap}{\code{signature(object =
	"hypredGenome")}: Modify the genetic map of the object}
    \item{hypredNewQTL}{\code{signature(object =
	"hypredGenome")}: Assign QTLs to the object.}
    \item{hypredRecombine}{\code{signature(object = "hypredGenome")}:
      Recombine two haploid genomes.}
    \item{hypredSNPdist}{\code{signature(object = "hypredGenome")}:
      Determine the distance and recombination frequency between two
      loci on a chromosome.}
    \item{hypredTruePerformance}{\code{signature(object =
	"hypredGenome")}: Determine the genotypic performance of
      individuals}
    \item{hypredTpGenomeSpecific}{\code{signature(object =
	"hypredGenome")}: Determine the genome specific genotypic performance
      of individuals }
    \item{show}{\code{signature(object = "hypredGenome")}: Print some
      summary information about the genome}
    \item{summary}{\code{signature(object = "hypredGenome")}: Print some
      summary information about the genome.}
	 }
       }

\author{
Frank Technow
}
\note{

  The number of SNP, QTL with additive effect, QTL with dominance
  effect, and perfect markers must be the same for each chromosome. This
  restriction helps allot in keeping track of the loci.
  
  The ID of a loci is its running number, where the first locus (the one with
  lowest map position) on chromosome 1 has ID = 1, the second locus (the
  one with the second lowest map position) on chromosome 1 has ID = 2
  and so on. When there are 10 loci per chromosome, then ID = 11 belongs
  to the first locus on chromosome 2.  
}

\section{Warning }{To assure a valid object, this object should only be created by the
  function \code{hypredGenome}, and modified with the functions
  \code{hypredNewMap} and \code{hypredNewQTL}.

} 


\seealso{ The function \code{\link{hypredGenome}} which creates the
  object, the function \code{\link{hypredNewMap}} which allows to modify
  the genetic map and the function \code{\link{hypredNewQTL}} which allows
  to assign QTL or modify already existing ones.}


\examples{
showClass("hypredGenome")
}
\keyword{classes}
