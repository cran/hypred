hypredDInv <- function(d, ## vector with genetic distances in Morgan
                 FUN) ## inverse mapping function to use ("Kosambi","Haldane")
  {
    if(FUN == "Kosambi")
      {
        r <- 0.5 * (exp(2*d) - exp(-2*d))/(exp(2*d) + exp(-2*d))
      }
    else if(FUN == "Haldane")## use Haldane
      {
        r <- 0.5 * (1 - exp(-2*d))
      }
    return(r)
  }
