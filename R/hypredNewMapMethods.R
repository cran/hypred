setMethod("hypredNewMap",
          signature(object = "hypredGenome"),
          function(object, new.map)
          {

            ##-------------------------------------------------##
            ## validity checks

            ## no negative elements

            if(any(new.map < 0))
              {
                stop("There must be no negative map positions!\n")
              }
            
            ## number of total loci must be identical to the original object
            
            if(!identical(object@num.chr * (object@num.snp.chr + object@num.add.qtl.chr),
                          length(new.map)))
              {
                stop("The length of new.map must match the number of loci in the object!\n")
              }


            for(i in 1:object@num.chr)
              {
                ## extract the map of the chromosome
                chrom.map <- new.map[object@chr.break.pts[i,1]:object@chr.break.pts[i,2]]

                ## positions within a chromosome must be in strictly
                ## increasing order!
                
                ## if the diff is zero or negative, the positions are not strictly increasing
                if(any(diff(chrom.map) <= 0))
                  {
                    stop("Map positions within a chromosome must be stricty increasing!\n")
                  }

                ## the maximum map position must not be beyond the
                ## length of the chromosome given in the object

                if(max(chrom.map) > object@len.chr[i])
                  {
                    stop(paste("The map positions in chromosome ", i,
                               " exceed its length!\n"))
                  }
              }
                ## everything ok?, than copy and modify the object and
                ## return it

                object.mod <- object

                ## change the map
                object.mod@pos.snp <- new.map

                ## if there are qtl, their positions must be changed
                ## as well

                if(object@num.add.qtl.chr > 0)
                object.mod@pos.add.qtl$M <- new.map[object@pos.add.qtl$ID]

                if(object@num.dom.qtl.chr > 0)
                object.mod@pos.dom.qtl$M <- new.map[object@pos.dom.qtl$ID]


                ## if there are blocks, give a warning that these
                ## positions will not be changed and that the slots of
                ## block.info will be set to NULL
                if(identical(is.null(object@block.info$num.blocks.chr),FALSE))
                  {
                    warning(paste("\nModifying the positions of SNP blocks is not supported!\n",
                                  "All slots storing block information are set to NULL!\n"))
                    
                    object.mod@block.info <- list(center = NULL,
                                                  end = NULL,
                                                  num.blocks.chr = NULL)


                    object.mod@pos.snp.block <- numeric()
                    
                  }
                
                return(object.mod)
                
          }
          )

