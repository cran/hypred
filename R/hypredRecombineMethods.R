setMethod("hypredRecombine",
          signature(object = "hypredGenome"),
          function(object, genomeA, genomeB, mutate, mutation.rate.snp, mutation.rate.qtl, block)
          {
            ## number of chromosomes
            num.chr <- object@num.chr
            
            ## number of SNP + addQTL per chromosome = total number of SNP
            
            tot.num.SNP <- object@num.snp.chr + object@num.add.qtl.chr

            if(identical(block, FALSE))
              {
                map <- object@pos.snp
              }
            else
              {
                if(is.null(object@block.info$num.blocks.chr))
                  {
                    stop("block = TRUE not possible because num.blocks.chr was NULL in hypredGenome")
                  }
                map <- object@pos.snp.block
              }
            
            ## validity checks ##

            if(length(genomeA) != tot.num.SNP * num.chr)
              {
                stop(paste("\nThe length of the genomeA input vector doesn't comply with the\n",
                           "number of SNP and QTL in the hypredGenome object\n"))
              }

            if(length(genomeB) != tot.num.SNP * num.chr)
              {
                stop(paste("\nThe length of the genomeB input vector doesn't comply with the\n",
                           "number of SNP and QTL in the hypredGenome object\n"))
              }

                
            ## vector that will store the recombined haploid

            hapVec <- integer(tot.num.SNP * num.chr)

            hapVec <-.C(meiosisFUNallChr,
                        as.integer(genomeA),
                        as.integer(genomeB),
                        newHaploid = integer(tot.num.SNP * num.chr),
                        as.double(map),
                        as.integer(num.chr),
                        as.double(object@len.chr),
                        as.integer(tot.num.SNP)
                        )$newHaploid   

            ## mutation part
            if(identical(mutate, TRUE))
              {
                if((mutation.rate.snp < 0) | (mutation.rate.qtl < 0))
                  stop("\nThe mutation rates must not be negative!\n")

                ##---------------------------------------##
                ## mutation of SNPs
                
                ## number of mutations in meiotic product follows a
                ## poison distribution with rate mutation.rate *
                ## number of SNP 
                num.mut <- rpois(1,mutation.rate.snp * object@num.snp.chr * num.chr)
                
                if(num.mut > 0) ## if there were mutations
                  {

                    ## if there are more mutations than SNP, allow
                    ## back mutations, by setting the argument replace
                    ## to TRUE in function sample. However, this
                    ## should happen only when the number of SNPs is
                    ## extremely small (i.e. 1,2,3)
                    back.mutate <- FALSE
                    
                    if(num.mut > (object@num.snp.chr * object@num.chr))
                      {
                        back.mutate <- TRUE
                      }
                      
                    ## - locations of the mutation are uniformely distributed
                    if(object@num.add.qtl.chr > 0)
                      {
                        loc.mut <- sample((1:(tot.num.SNP * num.chr))[-object@pos.add.qtl$ID],
                                          num.mut,
                                          replace = back.mutate)
                      }
                    else
                      {
                        loc.mut <- sample(1:(tot.num.SNP * num.chr),
                                          num.mut,
                                          replace = back.mutate)
                      }
                    
                    ## change 1 to 0 or 0 to 1
                    hapVec[loc.mut] <- (hapVec[loc.mut] - 1)^2    
                  }

                ##------------------------------------------##
                ## mutation of QTL

                ## number of mutations in meiotic product follows a
                ## poison distribution with rate mutation.rate.qtl *
                ## number of QTL
                num.mut <- rpois(1,mutation.rate.qtl * object@num.add.qtl.chr * num.chr)
                
                if(num.mut > 0) ## if there were mutations
                  {

                    ## if there are more mutations than QTL, allow
                    ## back mutations, by setting the argument replace
                    ## to TRUE in function sample. 
                    back.mutate <- FALSE
                    
                    if(num.mut > (object@num.add.qtl.chr * object@num.chr))
                      {
                        back.mutate <- TRUE
                      }
                      
                    ## - locations of the mutation are uniformely distributed
                    if(object@num.add.qtl.chr > 0)
                      {
                        ## from examples of ?sample
                        resample <- function(x, ...) x[sample.int(length(x), ...)]

                        loc.mut <- resample(object@pos.add.qtl$ID,
                                            num.mut,
                                            replace = back.mutate)

                        ## change 1 to 0 or 0 to 1
                        hapVec[loc.mut] <- (hapVec[loc.mut] - 1)^2    
                      }
                  } ## end if num.mut > 0

              } ## end mutation part

            
            hapVec <- matrix(hapVec,nrow = 1)
            colnames(hapVec) <- 1:(tot.num.SNP * num.chr)
            
            return(hapVec)
          }
          )



