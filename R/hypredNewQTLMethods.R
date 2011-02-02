setMethod("hypredNewQTL",
          signature(object = "hypredGenome"),
          function(object,
                   new.id.add = NULL,
                   new.id.dom = NULL,
                   new.id.per.mar = NULL,
                   new.eff.add = NULL,
                   new.eff.dom = NULL)
          {
            
            ## number of loci in object
            num.loci <- object@num.chr * (object@num.snp.chr + object@num.add.qtl.chr)

            
            ##-------------------------------------------------##
            ## validity checks
            ##-------------------------------------------------##

            
            ## there need to be additive QTL
            if(is.null(new.id.add))
              {
                stop("There needs to be at least one QTL\n")  
              }
            
                        
            ## ensure correct types

            new.id.add <- as.integer(new.id.add)
            new.eff.add <- as.numeric(new.eff.add)

            
            if(!is.null(new.id.dom))
              {    
                new.id.dom <- as.integer(new.id.dom)
              }
            if(!is.null(new.id.per.mar))
              {    
                new.id.per.mar <- as.integer(new.id.per.mar)
              }

            if(!is.null(new.eff.dom))
              {    
                new.eff.dom <- as.numeric(new.eff.dom)
              }

            ## ids of the additive QTL must be within the id range of
            ## SNP in the object
            
            if(!all(new.id.add %in% 1:num.loci))
              {
                stop("Not all QTL ids are within the id range of the SNP in the object!\n")
              }

            ## IDs must be unique
            if(!identical(length(new.id.add), length(unique(new.id.add))))
              {
                stop("The IDs of the qtl must be unique!\n")
              }

            ## IDs must be in strictly increasing order
            if(any(diff(new.id.add) < 0))
              {
                stop("The IDs of the qtl must be in strictly increasing order\n")
              }

            
            
            ##--- dominance qtl ---##


            if(!is.null(new.id.per.mar))
              {
                ## ids of the dominance QTL must be within the id range of
                ## newly defined additive QTL, i.e. only for QTL with
                ## additive effects may dominance effects be defined
                if(!all(new.id.dom %in% new.id.add))
                  {
                    stop(paste("Not all QTL with dominance effect have an id that\n",
                               "matches the ids of the QTL with additive effect!\n"))
                  }

                ## each ID can only appear once
                if(!identical(length(new.id.dom), length(unique(new.id.dom))))
                  {
                    stop("IDs of dominance qtl must be unique!\n")
                  }

                ## IDs must be in strictly increasing order
                if(any(diff(new.id.dom) < 0))
                   {
                     stop("The IDs of the dominance qtl must be in strictly increasing order\n")
                   }
              }

            ##--- perfect markers ---##

            if(!is.null(new.id.per.mar))
              {
                ## ids of the perfect marker must match the ids of the QTL 
                if(!all(new.id.per.mar %in% new.id.add))
                  {
                    stop(paste("ids of the perfect marker must match the ids of the QTL\n"))
                  }
                    
                ## each ID can only appear once
                if(!identical(length(new.id.per.mar), length(unique(new.id.per.mar))))
                  {
                    stop("Each ID of an perfect marker can only appear once!\n")
                  }
                
                ## IDs must be in strictly increasing order
                if(any(diff(new.id.per.mar) < 0))
                   {
                     stop("The IDs of the perfect markers must be in strictly increasing order\n")
                   }
              }
                
            ## lengths of vectors with effects must match the lengths
            ## of the vectors with ids

            if(!identical(length(new.id.add), length(new.eff.add)))
              {
                stop(paste("The number of additive effects doesn't match the number of QTL\n"))
              }

            if(!identical(length(new.id.dom), length(new.eff.dom)))
              {
                stop(paste("The number of dominance effects doesn't match the number of QTL with dominance effects\n"))
              }

                    
            ##--- this next part is to ensure that the numbers of     ---##
            ##--- purely additive qtl, qtl with dominance and perfect ---##
            ##--- markers are the same for all chromosomes            ---##

            ##--- This proceeds by calculating the expected number of  ---##
            ##--- qtl per chromosome and then checking if the new ids  ---##
            ##--- for qtl or marker on a given chromosome (extracted   ---##
            ##--- from the id vectors according to the expected number ---##
            ##--- of qtl/markers on chromosome) match the ids of snp   ---##
            ##--- on that chromosome                                   ---##
            
            ## number of qtl must be a multiple of the number of chromosomes
            if(!((length(new.id.add) %% object@num.chr ) == 0))
              {
                stop("The number of qtl must be a multiple of the number of chromosomes!\n")
              }
            
            ## calculate the number of additive qtl per chromosome
            num.add.qtl.chr.new <- as.integer(length(new.id.add)/object@num.chr)


            ## check if there are exactly num.add.qtl.chr.new per chromosome

            ## index variables
            first.locus.chr <- 1 # ID of first locus on the chromosome
            first.qtl.chr <- 1   # ID of first qtl on the chromosome

            num.loci.chr <- num.loci / object@num.chr ## number of loci on chromosome

            ## loop over chromosomes
            for(i in 1:object@num.chr)
              {
                ## IDs of the newly assigned qtl for this chromosome
                new.qtl.chr <- new.id.add[first.qtl.chr : (first.qtl.chr + num.add.qtl.chr.new - 1)]
                ## IDs of all the loci on the chromosome
                loci.chr <- first.locus.chr : (first.locus.chr + num.loci.chr - 1)

                ## if not all the QTLs have an ID that matches a ID of
                ## a SNP on the chromosome
                
                if(!all(new.qtl.chr %in% loci.chr))
                  {
                    stop("The number of qtl must be the same for all chromosomes!\n")
                  }

                ## increment the index variables
                first.locus.chr <- first.locus.chr + num.loci.chr 
                first.qtl.chr <- first.qtl.chr + num.add.qtl.chr.new
    
              }


            ##-- if there are qtl with dominance, do the same for them ---##

            if(!is.null(new.id.dom))
              {

                ## number of qtl must be a multiple of the number of chromosomes
                if(!((length(new.id.dom) %% object@num.chr ) == 0))
                  {
                    stop("The number of qtl with dominance must be a multiple of the number of chromosomes!\n")
                  }

                ## calculate the number of additive qtl per chromosome
                num.dom.qtl.chr.new <- as.integer(length(new.id.dom)/object@num.chr)


                ## check if there are exactly num.dom.qtl.chr.new per chromosome

                ## index variables
                first.locus.chr <- 1 # ID of first locus on the chromosome
                first.qtl.chr <- 1 # ID of first qtl on the chromosome

                num.loci.chr <- num.loci / object@num.chr ## number of loci on chromosome

                ## loop over chromosomes
                for(i in 1:object@num.chr)
                  {
                    ## IDs of the newly assigned qtl for this chromosome
                    new.qtl.chr <- new.id.dom[first.qtl.chr : (first.qtl.chr + num.dom.qtl.chr.new - 1)]
                    ## IDs of all the loci on the chromosome
                    loci.chr <- first.locus.chr : (first.locus.chr + num.loci.chr - 1)

                    ## if not all the QTLs have an ID that matches a ID of
                    ## a SNP on the chromosome
                
                    if(!all(new.qtl.chr %in% loci.chr))
                      {
                        stop("The number qtl with dominance must be the same for all chromosomes!\n")
                      }

                    ## increment the index variables
                    first.locus.chr <- first.locus.chr + num.loci.chr 
                    first.qtl.chr <- first.qtl.chr + num.dom.qtl.chr.new
                  }
              } ## end if dom

            
            ##--- if there are perfect markers, do the same for them ---##
                
            if(!is.null(new.id.per.mar))
              {
                ## number of perfect markers must be a multiple of the
                ## number of chromosomes
                if(!((length(new.id.per.mar) %% object@num.chr ) == 0))
                  {
                    stop("The number of perfect markers must be a multiple of the number of chromosomes!\n")
                  }

                ## calculate the number of perfect markers per chromosome
                num.per.mar.chr.new <- as.integer(length(new.id.per.mar)/object@num.chr)

                ## check if there are exactly num.per.mar.chr.new per chromosome

                ## index variables
                first.locus.chr <- 1 # ID of first locus on the chromosome
                first.per.mar.chr <- 1 # ID of first perfect marker on the chromosome

                num.loci.chr <- num.loci / object@num.chr ## number of loci on chromosome

                ## loop over chromosomes
                for(i in 1:object@num.chr)
                  {
                    ## IDs of the newly assigned perfect markers for this chromosome
                    new.per.mar.chr <- new.id.per.mar[first.per.mar.chr : (first.per.mar.chr + num.per.mar.chr.new - 1)]
                    ## IDs of all the loci on the chromosome
                    loci.chr <- first.locus.chr : (first.locus.chr + num.loci.chr - 1)

                    ## if not all the perfect markers have an ID that matches a ID of
                    ## a SNP on the chromosome
                
                    if(!all(new.per.mar.chr %in% loci.chr))
                      {
                        stop("The number of perfect markers must be the same for all chromosomes!\n")
                      }

                    ## increment the index variables
                    first.locus.chr <- first.locus.chr + num.loci.chr 
                    first.per.mar.chr <- first.per.mar.chr + num.per.mar.chr.new
    
                  } ## end if per.mar
              }     ## end for

                    
            ##--- end of this part ---##

            
            
            ## everything ok?, than copy and modify the object and
            ## return it

            object.mod <- object

            ##--- QTL with additive effect ---##
               
            ## change num.add.qtl.chr
            object.mod@num.add.qtl.chr <- num.add.qtl.chr.new
            
            ## change the ids of the additve QTL
            object.mod@pos.add.qtl$ID <- new.id.add

            ## change the map positions of the additve QTL
            object.mod@pos.add.qtl$M <- object@pos.snp[new.id.add]

            ## change the additive effects
            object.mod@add.and.dom.eff <- list(add = new.eff.add,
                                               dom = NULL)
               

            ##--- QTL with dominance ---##

               
            ## if new.id.dom is not NULL
            if(!is.null(new.id.dom))
              {
                ## change num.per.mar.chr
                object.mod@num.dom.qtl.chr <- num.dom.qtl.chr.new
                 
                ## set pos.dom.qtl 
                object.mod@pos.dom.qtl$ID <- new.id.dom
                object.mod@pos.dom.qtl$M  <- object@pos.snp[new.id.dom]
                 
                ## set the dominance effects
                object.mod@add.and.dom.eff$dom <- new.eff.dom
              }
            else
              {
                ## set num.dom.qtl.chr in the object to zero
                object.mod@num.dom.qtl.chr <- as.integer(0)
                 
                ## set pos.dom.qtl in the object to NULL
                object.mod@pos.dom.qtl <- list(ID = NULL,
                                               M = NULL)
                 
                ## dominance effects were set to NULL already
                
              }


               
            ##--- perfect markers ---##

            ## if new.id.per.mar is not NULL
            if(!is.null(new.id.per.mar))
              {
                ## change num.per.mar.chr
                object.mod@num.per.mar.chr <- num.per.mar.chr.new
                 
                ## change the IDs 
                object.mod@per.mar.id <- new.id.per.mar
              }
            else
              {
                ## set num.per.mar.chr in the object to zero
                object.mod@num.per.mar.chr <- as.integer(0)
                 
                ## set per.mar.id in the object to NULL
                object.mod@per.mar.id <- as.integer(0)
              }


            ##--- The number of loci always remains constant.          ---##
            ##--- Hence, change num.snp.chr according to the number of ---##
            ##--- qtl per chr                                          ---##

            object.mod@num.snp.chr <- as.integer((num.loci - length(new.id.add)) / object@num.chr)

            
            return(object.mod)
          }
          )


