setMethod("hypredCode",
          signature(object = "hypredGenome"),
          function(object, genotypes, DH, type)
          {

            ## validity checks

            ## number of snp marker
            num.snp <- (object@num.snp.chr + object@num.per.mar.chr) * object@num.chr
            
            if(!identical(ncol(genotypes),
                          as.integer((object@num.snp.chr + object@num.add.qtl.chr) *
                                     object@num.chr)))
              {
                stop("The number of columns in genotypes doesn't agree with the total number of loci")
              }
            
            ## give genotypes colum names             
            colnames(genotypes) <- as.character(1:as.integer((object@num.snp.chr +
                                                              object@num.add.qtl.chr) *
                                                             object@num.chr))
            
            ## remove the QTL-SNPs, exept those that correspond to
            ## perfect markers

            if(object@num.per.mar.chr < object@num.add.qtl.chr)
              {
                genotypes <-
                  genotypes[,-object@pos.add.qtl$ID[!(object@pos.add.qtl$ID %in% object@per.mar.id)]]
              }
            else 
              {
                ## leave genotypes unchanged, since all QTL correspond
                ## to perfect markers (because num.add.qtl.chr =
                ## num.per.mar.chr) The case that num.per.mar.chr >
                ## num.add.qtl.chr is ruled out by the function
                ## hypredGenome
              }
            
            ## DH = FALSE
            

            if(identical(DH, FALSE))
              {
                N <- nrow(genotypes)
                
                if(!identical(as.integer(N %% 2), as.integer(0)))
                  {
                    stop(paste("When DH = FALSE, the number of rows in the actual\n",
                               "argument to genotype must be a even number!"))
                  }
                
                N <- N/2
              }

            ## Xu2003
            if(identical(DH, FALSE) & identical(type,"Xu2003"))
              {
               
                res <- .C(hypredCode_FUN_dom,
                          as.integer(genotypes),
                          X = double(N * num.snp),
                          W = double(N * num.snp),
                          as.integer(num.snp),
                          as.integer(N)
                          )

                X <- matrix(res$X, nrow = N, ncol = num.snp)
                ## if colnames(genotypes) = NULL, colnames(X) = NULL as well
                colnames(X) <- colnames(genotypes)
                W <- matrix(res$W, nrow = N, ncol = num.snp)
                colnames(W) <- colnames(genotypes)
                return(cbind(X,W))
                     }
            
            else  if(identical(DH, FALSE) & identical(type,"012"))
              {
                first <- 1
                X <- matrix(nrow = N, ncol = num.snp)
                colnames(X) <- colnames(genotypes)
                for(i in 1:N)
                  {
                    X[i,] <- colSums(genotypes[c(first, first + 1),])
                    first <- first + 2 ## first homologe of individual i+1
                  }
                
                return(X)
              }

            else if(identical(DH, FALSE) & identical(type,"-101"))
              {
                first <- 1
                X <- matrix(nrow = N, ncol = num.snp)
                colnames(X) <- colnames(genotypes)
                for(i in 1:N)
                  {
                    X[i,] <- colSums(genotypes[c(first, first + 1),])
                    first <- first + 2 ## first homologe of individual i+1
                  }
                
                return(X-1)
              }
            
            else if(identical(DH, FALSE) & identical(type,"genome.specific"))
              {
                ## [X1,X2]
                X1X2 <- matrix(t(genotypes), nrow = N, byrow = TRUE)
                colnames(X1X2) <- rep(colnames(genotypes),2)
                ## W1W2
                
                W1W2 <- codeW.genomeSpecific(genotypes, N, num.snp)
                colnames(W1W2) <- rep(colnames(genotypes),2)
                return(cbind(X1X2,W1W2))
              }

            
            ## if DH = TRUE
            
            if(identical(DH, TRUE) & identical(type,"012"))
              {
                return(genotypes * 2)
              }

            if(identical(DH, TRUE) & identical(type,"-101"))
              {
                return(genotypes * 2 - 1)
              }

            if(identical(DH, TRUE) & identical(type,"Xu2003"))
              {
                stop("\nUsing a DH population to estimate dominance effects is not a good idea!\n")
              }

            if(identical(DH, TRUE) & identical(type,"genome.specific"))
              {
                stop("\nThis can't be done for DH!\n")
              }
            
        } ## end function
)
