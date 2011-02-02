setMethod("hypredTpGenomeSpecific",
          signature(object = "hypredGenome"),
          function(object,
                   genotypes,
                   add.eff.G1 = NULL,
                   add.eff.G2 = NULL,
                   dom.eff.G1 = NULL,
                   dom.eff.G2 = NULL
                   )
          {
            
            if(identical(object@num.add.qtl.chr, as.integer(0)))
              {
                stop("\nThere are no QTL specified in the object!\n")
              }
              
            if(is.null(add.eff.G1))
              {
                add.eff.G1 <- object@add.and.dom.eff$add
              }
            if(is.null(add.eff.G2))
              {
                add.eff.G2 <- object@add.and.dom.eff$add
              }
            if(is.null(dom.eff.G1))
              {
                dom.eff.G1 <- object@add.and.dom.eff$dom
              }
            if(is.null(dom.eff.G2))
              {
                dom.eff.G2 <- object@add.and.dom.eff$dom
              }

            ## validity checks

            if(!identical(ncol(genotypes),
                          as.integer((object@num.snp.chr + object@num.add.qtl.chr) *
                                     object@num.chr)))
              {
                stop("The number of columns in genotypes doesn't agree with the total number of loci")
              }


            if(!is.null(add.eff.G1) & (length(add.eff.G1) != (object@num.add.qtl.chr * object@num.chr)))
              {
                stop(paste("The lengths of the vectors with new add. effects must\n",
                           "be equal to the number of add. QTL in the object"))
              }
            
            if(!is.null(add.eff.G2) & (length(add.eff.G2) != (object@num.add.qtl.chr * object@num.chr)))
              {
                stop(paste("The lengths of the vectors with new add. effects must\n",
                           "be equal to the number of add. QTL in the object"))
              }


            if(identical(object@num.dom.qtl.chr, as.integer(0)))
              {
                if(!is.null(dom.eff.G1) | !is.null(dom.eff.G2))
                  
                  stop("\nThere were no dom. QTL specified in the object\n")

              }
              
            
            if(!is.null(dom.eff.G1) & (length(dom.eff.G1) != (object@num.dom.qtl.chr * object@num.chr)))
              {
                stop(paste("The lengths of the vectors with new dom. effects must\n",
                           "be equal to the number of dom. QTL in the object"))
              }

            if(!is.null(dom.eff.G2) & (length(dom.eff.G2) != (object@num.dom.qtl.chr * object@num.chr)))
              {
                stop(paste("The lengths of the vectors with new add. effects must\n",
                           "be equal to the number of dom. QTL in the object"))
              }




            ## extract the SNPs that are additive QTL 
            QTL.genotype.add <- genotypes[,object@pos.add.qtl$ID, drop = FALSE]

            ## find number of individuals

            N <- nrow(QTL.genotype.add)
              
            if(N %% 2 != 0){
              stop(paste("For hybrids, the number of rows in the actual\n",
                         "argument to genotype must be an even number!"))
            }

            N <- N/2


            ## additive effects y1(add) = X[1,]b1 + X[2,]b2
            ##
            ##where b1 and b2 are the effects for population 1 resp 2.

            ## sum of effects of G1
            true.add.G1 <- QTL.genotype.add[((1:(N*2)) %% 2) != 0,] %*% matrix(add.eff.G1, ncol = 1)
            ## sum of effects of G2
            true.add.G2 <- QTL.genotype.add[((1:(N*2)) %% 2) == 0,] %*% matrix(add.eff.G2, ncol = 1)

            ## sum
            true.performance.add <- true.add.G1 + true.add.G2
                                    
            true.performance.dom <- rep(0,N)
              
            if(!identical(object@num.dom.qtl.chr, as.integer(0))){
              

              ## extract the SNPs that have dominance effects
              ## QTL 
              QTL.genotype.dom <- genotypes[,object@pos.dom.qtl$ID, drop = FALSE]

              W1W2 <- codeW.genomeSpecific(QTL.genotype.dom, N,
                                           object@num.dom.qtl.chr * object@num.chr)
              
              ## W1W2 * d1d2
              true.performance.dom <- W1W2 %*% c(dom.eff.G1, dom.eff.G2)

            } ## end if num.dom.qtl.chr != 0

            true.performance <- as.vector(true.performance.add) +
              as.vector(true.performance.dom)
              

            dimnames(true.performance) <- c(NULL,NULL)

            return(as.matrix(true.performance))
            
          } ## end function
          )
