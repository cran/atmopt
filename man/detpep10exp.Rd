\name{detpep10exp}
\alias{detpep10exp}
\docType{package}
\title{
  DetPep10Exp discrete test function
}
\description{
  A discrete test function constructed from the modified exponential function in Dette and Pepelyshev (2010).
}
\usage{
  detpep10exp(xx,ntimes,nlev)
}
\arguments{
  \item{xx}{A \eqn{p}-vector for input factors.}
  \item{ntimes}{Number of duplications for the function (base function is 3D).}
  \item{nlev}{A \eqn{p}-vector corresponding to the number of levels for each factor(discretized on equally-spaced intervals).}
}
\examples{
xx <- c(1,2,1,2,1,2) #input factors
nlev <- rep(4,length(xx)) #number of levels for each factor
ntimes <- length(xx)/3 #base function is in 2D, so duplicate 2 times
detpep10exp(xx,ntimes,nlev)
}
