detpep10exp <- function(xx,ntimes,nlev)
{
  #transform to unit cube
  xx <- trans.uh(matrix(xx,nrow=1),nlev,levs=NULL)

  #compute function
  y <- 0.0
  for (ll in 1:ntimes){
    x1 <- xx[3*ll-2]
    x2 <- xx[3*ll-1]
    x3 <- xx[3*ll]

    term1 <- exp(-2/(x1^1.75))
    term2 <- exp(-2/(x2^1.5))
    term3 <- exp(-2/(x3^1.25))
    term4 <- 0.01*x1*x2*x3

    y <- y + 100 * (term1 + term2 + term3 + term4)
  }
  return(y)
}

camel6 <- function(xx,ntimes,nlev)
{

  #transform to unit cube
  xx <- trans.uh(matrix(xx,nrow=1),nlev,levs=NULL)

  #compute function
  runsum <- 0.0
  for (ll in 1:ntimes){
    x1 <- 4*xx[2*ll-1]-2
    x2 <- 2*xx[2*ll]-1

    term1 <- (4-2.1*x1^2+(x1^4)/3) * x1^2
    term2 <- x1*x2
    term3 <- (-4+4*x2^2) * x2^2

    runsum <- runsum + term1 + term2 + term3

  }
  return(runsum)
}
