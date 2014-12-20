\name{hypred-package}
\alias{hypred-package}
\alias{hypred}
\docType{package}
\title{
hypred: Simulation of Genomic Data in Applied Genetics
}

\description{ The package \pkg{hypred} is intended for simulating high
  density genomic data. The meiosis is simulated under a count-location
  model, the genetic structure can include additive and dominance
  effects. The low level nature of the package provides great
  flexibility in creating all kinds of populations.  }

\details{
\tabular{ll}{
Package: \tab hypred\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2014-12-20\cr
License: \tab GNU General Public License Version 2, June 1991\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}

The function \code{hypredGenome} is used to define the genome
parameters, number and length of chromosomes, and number of SNP loci. It
also simulates a genetic map. This functions returns an object that is
of class \code{hypredGenome}. The function \code{hypredNewMap} is used
to modify the genetic map, and the function \code{hypredNewQTL} to
assign QTL. The function \code{hypredFounder} will create two homozygous
founder lines according to the object returned by hypredGenome. The
meiosis is simulated by the function \code{hypredRecombine}. This
function takes two haploid chromosome sets (as found in a sex cell) and
returnes a recombined haploid chromosome set (i.e. the gamete). The
meiosis is simulated according to the count-location model, simulation
of mutations is possible as well. The functions
\code{hypredTruePerformance} and \code{hypredTpGenomeSpecific} return
the true genotypic values of individuals given their
genotypes. \code{hypredTruePerformance} uses the QTL effects found in
the \code{hypredGenome} object, \code{hypredTpGenomeSpecific} allows to specify
new and genome specific effects. Various types of design matrices (to be
used in modeling algorithms) are created from the genotypes by the
function \code{hypredCode}. All functions mentioned above are generic,
with methods defined for the \code{hypredGenome} object.}

\author{Frank Technow

  Maintainer: Frank Technow \email{ Frank.Technow@gmx.net }
}

\keyword{ package }
\keyword{datagen}

