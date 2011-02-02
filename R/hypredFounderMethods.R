setMethod("hypredFounder",
          signature(object = "hypredGenome"),
          function(object, prob.snp)
          {

            ## assure that 0 <= prob.snp <= 1

            if((prob.snp < 0) | (prob.snp > 1))
              {
                stop("prob.snp must be within [0, 1]!\n")
              }
            
            founder1 <- sample(c(1,0),
                               size = object@num.chr * (object@num.snp.chr + object@num.add.qtl.chr),
                               replace = TRUE,
                               prob = c(prob.snp, 1 - prob.snp)
                               )
            founder2 <- as.integer(!founder1)

            
            names(founder1) <- 1:(object@num.chr * (object@num.snp.chr + object@num.add.qtl.chr))
            
            return(rbind(founder1,
                         founder2))
          }
          )


