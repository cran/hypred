hypredGenome <- function(num.chr,
                         len.chr,
                         num.snp.chr,
                         num.blocks.chr = NULL
                         )
  {

    ## assure correct types

    num.chr <- as.integer(num.chr)
    len.chr <- as.numeric(len.chr)
    num.snp.chr <- as.integer(num.snp.chr)
    
    if(!is.null(num.blocks.chr))
      {
        num.blocks.chr <- as.integer(num.blocks.chr)
      }

    ##--------- validity checks ------------##

    ## no negative arguments

    if(any(c(num.chr, len.chr, num.snp.chr, num.blocks.chr) < 0))
      {
        stop("\nThere must be no negative arguments!\n")
      }
    
    ## length of len.chr must comply with num.chr

    if(!identical(num.chr, length(len.chr)))
       {
         stop(paste("The number of chromosomes (num.chr) and the length\n",
                    "of the vector of chromosome lenghts (len.chr) must be equal!",
                    sep = ""
                    ))
       }


    ## if (num.blocks.chr/num of SNP) < 0.1, stop

    if(!is.null(num.blocks.chr))
      {
        if((num.blocks.chr/num.snp.chr) > 0.100001)
          {
            stop(paste("\nWhen non NULL, the argument num.blocks.chr must such that,\n",
                       "(num.blocks.chr/num.snp.chr) <= 0.1", sep = ""))
          }
      }
    
    ## create matrix that stores the positions of the SNP

    pos.mat <- matrix(nrow = num.chr,
                     ncol = num.snp.chr)


    ## sample the map positions from a uniform distribution [0,len.chr]
    for(i in 1:num.chr)
      {
        pos.mat[i,] <- sort(runif(num.snp.chr,
                                  min = 0,
                                  max = len.chr[i]))
      }

    ## turn pos.mat into vector

    pos.vec <- as.vector(t(pos.mat))

    ## postitions of blocks

    block.pos.vec <- numeric()
    center.SNP <- NULL
    end.SNP <- NULL

    if(!is.null(num.blocks.chr)) ## only if non NULL
      {
    ## same as above, but with positions of SNP blocks instead of
    ## positions of individual SNPs

    ## the position of a SNP block is the position of its center SNP
    ## This center is found as:
    ##
    ## block1: floor(num of SNP / 2*num of blocks)
    ## block2: floor(num of SNP / 2*num of blocks) + 1*floor(num of SNP / num of blocks)
    ## block3: floor(num of SNP / 2*num of blocks) + 2*floor(num of SNP / num of blocks)
    ## ...

    ## find center SNP
    
    offset <- floor(num.snp.chr/(num.blocks.chr * 2))
    increment <- floor(num.snp.chr/num.blocks.chr)

    center.SNP <- integer(num.blocks.chr)

    for(i in 1:num.blocks.chr)
      {
        center.SNP[i] <- offset + increment * (i - 1)
      }

    ## find the last SNP of a block
    
    end.SNP <- integer(num.blocks.chr)

    for(i in 1:(num.blocks.chr - 1))
      {
        end.SNP[i] <-  increment * i
      }
    
    end.SNP[num.blocks.chr] <- num.snp.chr ## last one just end of
                                           ## chromosome
                                          
    ## calculate the block positions

    block.pos.mat <- matrix(nrow = num.chr,
                            ncol = num.snp.chr)

    for(i in 1:num.chr)
      {
        start.SNP <- 1
        for(j in 1:num.blocks.chr)
          {
            block.pos.mat[i, start.SNP : end.SNP[j]] <-
              pos.mat[i,center.SNP[j]]
            
            start.SNP <- end.SNP[j] + 1
          }
        
      }

    ## turn block.pos.mat into vector

    block.pos.vec <- as.vector(t(block.pos.mat))
  }
    
    ## create matrix that stores the start and end points of a
    ## chromosome in the vectors
    
    ## for loop counter for the current chromosome
    current.chr <- 0
    chr.break.pts <- matrix(nrow = num.chr, ncol = 2)
    
    for(i in 1:num.chr)
      {
        ## first SNP of the current chromosome in the vectors
        chr.break.pts[i,1] <- current.chr * num.snp.chr + 1

        ## last SNP of the current chromosome in the vectors
        chr.break.pts[i,2] <- (current.chr + 1) * num.snp.chr
        
        current.chr <- current.chr + 1
      }


#### create and return an object of class hypredGenome
    
    res <- new("hypredGenome",
               num.chr = num.chr,
               len.chr = len.chr,
               num.snp.chr = num.snp.chr,
               pos.snp = pos.vec,
               pos.snp.block = block.pos.vec,
               block.info =
               list(center = center.SNP,
                    end = end.SNP,
                    num.blocks.chr = num.blocks.chr),
               chr.break.pts = chr.break.pts,
               num.add.qtl.chr = as.integer(0),
               num.dom.qtl.chr = as.integer(0),
               pos.add.qtl =
               list(ID = NULL,
                    M = NULL),
               pos.dom.qtl =
               list(ID = NULL,
                    M = NULL),
               add.and.dom.eff =
               list(add = NULL,
                    dom = NULL),
               num.per.mar.chr = as.integer(0),
               per.mar.id = as.integer(0))
  }
