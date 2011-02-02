setMethod("hypredSNPdist",
          signature(object = "hypredGenome"),
          function(object, chromosome, SNP1, SNP2, block = FALSE)
          {

            ## validy checks
            if(chromosome > object@num.chr)
              {
                stop(paste("Number of chromosomes given (",
                          chromosome,
                          ") \nis greater than the number of chromosomes defined in the object (",
                          object@num.chr, ")"))
              }

            if(SNP1 >= SNP2)
              {
                stop("SNP1 must have a smaller number than SNP2")
              }


            if(SNP2 > (object@num.snp.chr + object@num.add.qtl.chr))
              {
                stop(paste("Number of SNP2 (",
                           SNP2,
                           ") \nis greater than the number of SNP and QTL defined in the object (",
                           (object@num.snp.chr + object@num.add.qtl.chr), ")"))
              }
            
            
            if(identical(block, FALSE))
              {
                chrom <- object@pos.snp[object@chr.break.pts[chromosome,1] :
                                        object@chr.break.pts[chromosome,2] ]
              }
            else
              {
                if(is.null(object@block.info$num.blocks.chr))
                  {
                    stop("block = TRUE not possible because num.blocks.chr was NULL in hypredGenome")
                  }
                
                chrom <- object@pos.snp.block[object@chr.break.pts[chromosome,1] :
                                              object@chr.break.pts[chromosome,2] ]
              }
            
            dist <- diff(chrom[c(SNP1, SNP2)])
            rec <- hypredDInv(dist,"Haldane")
            
            return(rbind("d" = dist,
                         "r" = rec))
          }
          )
