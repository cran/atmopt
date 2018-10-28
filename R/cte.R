cte <- function(xx,alpha){
  #Conditional tail expectation
  return( mean(xx[xx<(quantile(xx,alpha)+1e-10)]) )
}

cte.fnl <- function(alphas){
  #Returns CTE functional list for each factor
  nfactor <- length(alphas)
  fn.lst <- vector("list",nfactor)
  for (j in 1:nfactor){
    fn.lst[[j]]$fn <- function(xx){cte(xx,alphas[j])}
    fn.lst[[j]]$alpha <- alphas[j]
  }
  return( fn.lst )
}
