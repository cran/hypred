setMethod("show",
          signature(object = "hypredGenome"),
          function(object)
          {
            sumRes <- list(num.chr = object@num.chr,
                           len.chr = object@len.chr,
                           num.snp.chr = object@num.snp.chr,
                           num.add.qtl.chr = object@num.add.qtl.chr,
                           num.dom.qtl.chr = object@num.dom.qtl.chr,
                           num.per.mar.chr = object@num.per.mar.chr,
                           sumAddEff = (object@add.and.dom.eff$add),
                           sumDomEff = (object@add.and.dom.eff$dom))

            ##
            cat("\nThe defined genome is characterized as follows:\n")
            cat(paste("\n",sumRes$num.chr," chromosomes of lengths ", sep = "") ,sumRes$len.chr," M\n", sep = " ")
            cat(sumRes$num.snp.chr," SNP markers plus ",sumRes$num.per.mar.chr," perfect SNP markers per chrom.\n", sep = "")
            cat(sumRes$num.add.qtl.chr," QTL with additive effect per chrom.\n", sep = "")
            cat("Of these, ",sumRes$num.dom.qtl.chr," also have (has) a dominance effect\n\n", sep = "")
            if(object@num.add.qtl.chr > 0)
              {
                cat("Summary of distribution of additive effects:\n\n")
                print(summary(sumRes$sumAddEff))
                cat("\n\n")
              }
            if(object@num.dom.qtl.chr > 0)
              {
                cat("Summary of distribution of dominance effects:\n\n")
                print(summary(sumRes$sumDomEff))
                cat("\n")
              }

          }
)
