setGeneric("summary")
setMethod("summary",
          signature(object = "hypredGenome"),
          function(object,showSummary = TRUE)
          {
            sumRes <- list(num.chr = object@num.chr,
                           len.chr = object@len.chr,
                           num.snp.chr = object@num.snp.chr,
                           num.add.qtl.chr = object@num.add.qtl.chr,
                           num.dom.qtl.chr = object@num.dom.qtl.chr,
                           num.per.mar.chr = object@num.per.mar.chr,
                           sumAddEff = (object@add.and.dom.eff$add),
                           sumDomEff = (object@add.and.dom.eff$dom),
                           chr.break.pts = object@chr.break.pts)
            ##
            if(identical(showSummary, TRUE))
              {
                show(object)
              }
            return(invisible(sumRes))
          }
          )
