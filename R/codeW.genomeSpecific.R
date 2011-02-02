codeW.genomeSpecific <- function(x, N, num.snp)
  {

    ## this turns all 1 in every first row to -1
    x <- x * rep(c(-1,1), N)

    ## sum the rows corresponding to one individual and
    ## (this will turn a (1,1) -> 0; (0,0) -> 0; (1,0) -> -1; and
    ## (0,1) -> 1

    first <- 1 ## the first homologe of individual i
    W <- matrix(nrow = N,
                ncol = num.snp)

    if(num.snp > 1){
      for(i in 1:N)
        {
          W[i,] <- colSums(x[c(first, first + 1),])
          first <- first + 2 ## first homologe of individual i+1
        }
    }
    else{
      for(i in 1:N)
        {
          W[i,] <- sum(x[c(first, first + 1),])
          first <- first + 2 ## first homologe of individual i+1
        }
    }

    W1W2 <- cbind(W,W)

    ## process W1: set all elements of W1 to zero that are 1,
    cols.of.W1 <- 1 : num.snp
              
    W1W2[,cols.of.W1][W1W2[,cols.of.W1] == 1] <- 0

    ## and square the result
    W1W2[,cols.of.W1] <- W1W2[,cols.of.W1]^2
    ## now an element of W1 is 1 only if the individual was
    ## 1 0 for this locus

    ## process W2: set all elements of W2 to zero that are -1

    cols.of.W2 <- (num.snp + 1) :
      (2 * num.snp)
              
    W1W2[,cols.of.W2][W1W2[,cols.of.W2] == -1] <- 0

    ## squaring is not necessary
              
    ## now an element of W2 is 1 only if the individual was
    ## 0 1 for this locus

    return(W1W2)

  }
