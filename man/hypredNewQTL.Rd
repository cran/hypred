\name{hypredNewQTL}
\docType{methods}
% Generic
\alias{hypredNewQTL}
% Methods
\alias{hypredNewQTL,hypredGenome-method}

\title{
Assign new QTLs to a hypredGenome object
}

\description{ A generic function to assign new QTLs to a
\code{"\linkS4class{hypredGenome}"} object }

\usage{
hypredNewQTL(object, \dots)
\S4method{hypredNewQTL}{hypredGenome}(object, new.id.add = NULL,
  new.id.dom = NULL, new.id.per.mar = NULL, new.eff.add = NULL,
  new.eff.dom = NULL)
}
\arguments{
  \item{object}{an object of a class that holds genomic information on map
    positions of loci, typically an \code{"\linkS4class{hypredGenome}"} object.}
  \item{new.id.add}{integer vector giving the IDs of the QTL with
    additive effects (hence all QTL).}
  \item{new.id.dom}{integer vector giving the IDs of the QTL defined
    with \code{new.id.add} that have also a dominance effect.}
  \item{new.id.per.mar}{integer vector giving the IDs of the QTL defined
    with \code{new.id.add} that are taken to be perfect markers.}
  \item{new.eff.add}{numeric vector giving the additive effects of the
    QTL defined with \code{new.id.add}.}
  \item{new.eff.dom}{numeric vector giving the dominance effects of the
    QTL defined with \code{new.id.dom}.}
  \item{\dots}{Methods may require further arguments.}
}
\value{
A \code{"\linkS4class{hypredGenome}"} object, modified accordingly.
}
\details{
The number of loci, as defined through the function
\code{hypredGenome} always remains constant. What can be changed is the
number of these that are QTL.

The number of QTL, as well as QTL with dominance effect, must be the
same for all chromosomes. Also the number of perfect markers must be the
same for all chromosomes. There are some validity checks included to
assure that. However, simply setting effects of some of the QTL to 0, is
an indirect way to have different numbers of QTL per chromosome.

Only QTL with an additive effect may have dominance effects
assigned. Hence, \code{new.id.dom} must be a subset of
\code{new.id.add}. However, setting the additive effect of a QTL with
dominance to zero is again an indirect way to have QTL with pure
dominance.  Of course, only QTL can be perfect markers, hence
\code{new.id.per.mar} must be a subset of \code{new.id.add} as well.


The vector given to \code{new.id.add} holds the IDs of the locis that
are assigned QTLs. The ID of a loci is its index in the vector of all
loci in the object. If, for example, the object hold 20 loci, then the
ID of the loci with lowest map postion of chromosome 1 is 1, the next
has ID 2 and so on. The ID of the loci with lowes map position of
chromosome 2 is 11. So if one wants to assign QTLs to the third locus on
chromosome 1 and the fifths on chromosome 2, the vector given to
\code{new.eff.add} needs to be \code{c(3, 15)}.

The same applies for the vectors given to \code{new.id.dom} and
\code{new.id.per.mar}.  

The effects are distributed to the QTL in the order in which they appear
in the vectors given to \code{new.eff.add} and \code{new.eff.dom}. So,
inline with above example, \code{new.eff.add = c(0.1, -0.1)} would
assign an additive effect of 0.1 to the QTL on chromosome 1 and -0.1 to the QTL
on chromosome 2.

If no dominance QTL or perfect markers are assigned, the corresponding
argument needs to be \code{NULL}.

The additive effects of the QTL are interpreted as the homozygous
effects (a) of the Falconer scale (i.e. half the difference between the
two homozygous genotypes). The dominance effects are interpreted as the
d of the Falconer scale (i.e. the difference between the heterozygous
genotype from the mean of the two homozygous genotypes).

The interpretation and treatment of QTL is similar to SNPs. It is
assumed that the QTL alleles are identical to each other, apart from one
causative SNP mutation. This mutation allows to identify them, so that a
1 codes for the mutated allele and a 0 for the unmutated allele. These
causative SNP are treated as unobserved by default, hence they are
excluded from the design matrices created by function
\code{hypredCode}. However, with the argument \code{new.id.per.mar}
causative SNP can be turned into observed, perfect markers that are
included in the design matrices.
}

\seealso{
The function \code{\link{hypredGenome}} which creates the original
object and the function \code{\link{hypredTruePerformance}} which
combines the QTL information into genotypic values of individuals.
}

\references{
  D. S. Falconer (1967)
  \emph{Introduction to Quantitative Genetics}. Oliver and Boyd.,
  London, 1967.
}


\examples{

## two chromosomes of length 1 M and 10 SNP (loci) per chromosome

no_QTL <- hypredGenome(2, c(1.0, 1.0), 10)

summary(no_QTL)

## assign a QTL to the third loci on chromosome 1 and the fifths on
## chromosome 2.

## The additive effect of the first QTL is 0.1, the additive effect of
## the second QTL is -0.1. The QTLs don't show dominance, both are
## perfect markers. 

with_QTL <- hypredNewQTL(no_QTL,
                         new.id.add = c(3, 15),
                         new.id.dom = NULL, ## default
                         new.id.per.mar = c(3, 15),
                         new.eff.add = c(0.1, -0.1),
                         new.eff.dom = NULL ## default
                        )

summary(with_QTL)

}

\keyword{methods}
\keyword{datagen}


