setMethod("hypredTruePerformance",
          signature(object = "hypredGenome"),
          function(object, genotypes, DH)
          {

            ## validity checks
            
            if(identical(object@num.add.qtl.chr, as.integer(0)))
              {
                stop("\nThere are no QTL specified in the object!\n")
              }

            if(!identical(ncol(genotypes),
                          as.integer((object@num.snp.chr + object@num.add.qtl.chr) *
                                     object@num.chr)))
              {
                stop("The number of columns in genotypes doesn't agree with the total number of loci")
              }

            
            ## extract the SNPs that are additive QTL (matrix X,
            ## X[i,j] is 1 if the jth SNP is present in the ith
            ## individual (if DH == TRUE) or in the ith homologe,
            ## where homologes (1,2) belong to individual 1, homologes
            ## (3,4) to individual 2 and so on.)
            QTL.genotype.add <- genotypes[,object@pos.add.qtl$ID, drop = FALSE]

            ## if DH, y = Xb * 2
            if(identical(DH,TRUE)){
              true.performance <- (QTL.genotype.add %*% object@add.and.dom.eff$add)*2
            }

            else{
              
              ## find number of individuals

              N <- nrow(QTL.genotype.add)
              
              if(N %% 2 != 0){
                stop(paste("When DH = FALSE, the number of rows in the actual\n",
                           "argument to genotypes must be a even number!"))
              }

              N <- N/2
              
              ## additive effects
              ## y1(add) = X[1,]b + X[2,]b
              
              ## indicator of individual
              fact <- rep(1:N, each = 2)
              ## sum of effects of one row (X[1,] resp. X[2,])
              true.add.one.row <- QTL.genotype.add %*% object@add.and.dom.eff$add
              ## X[1,]b + X[2,]b
              true.performance.add <- tapply(true.add.one.row,
                                             as.factor(fact),
                                             sum)


              true.performance.dom <- rep(0,N)
              
              if(!identical(object@num.dom.qtl.chr, as.integer(0))){
              
                ## dominance effects
                ## y1(dom) = W[1,]u

                ## extract the SNPs that have dominance effects
                ## QTL 
                QTL.genotype.dom <- genotypes[,object@pos.dom.qtl$ID, drop = FALSE]
              
                ## this turns all 1 in every second row to -1
                QTL.genotype.dom <- QTL.genotype.dom * rep(c(1,-1), N)

                ## sum the rows corresponding to one individual and
                ## square each sum (this will turn a (1,1) -> 0; (0,0)
                ## -> 0; (1,0) and (0,1) -> 1, hence W[i,j] will be 1 if
                ## the individual i is heterozygot for loci j and 0 if not.

                first <- 1 ## the first homologe of individual i
                W <- matrix(nrow = N,
                            ncol = object@num.dom.qtl.chr * object@num.chr)

                if((object@num.dom.qtl.chr * object@num.chr) > 1)
                  {
                    for(i in 1:N)
                      {
                        W[i,] <- colSums(QTL.genotype.dom[c(first, first + 1),])^2
                        first <- first + 2 ## first homologe of individual i+1
                      }
                  }
                else
                  {
                    for(i in 1:N)
                      {            
                        W[i,] <- sum(QTL.genotype.dom[c(first, first + 1),])^2
                        first <- first + 2 ## first homologe of individual i+1
                      }
                  }
                
                true.performance.dom <- W %*% object@add.and.dom.eff$dom

              } ## end if num.dom.qtl.chr != 0

              true.performance <- as.vector(true.performance.add) +
                as.vector(true.performance.dom)
              
            } ## end else !DH

            dimnames(true.performance) <- c(NULL,NULL)

            return(as.matrix(true.performance))
          } ## end function
          )

