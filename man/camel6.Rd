\name{camel6}
\alias{camel6}
\docType{package}
\title{
  Six-hump discrete test function
}
\description{
  A discrete test function constructed from the six-hump camel function in Ali et al. (2005).
}
\usage{
  camel6(xx,ntimes,nlev)
}
\arguments{
  \item{xx}{A \eqn{p}-vector for input factors.}
  \item{ntimes}{Number of duplications for the function (base function is 2D).}
  \item{nlev}{A \eqn{p}-vector corresponding to the number of levels for each factor(discretized on equally-spaced intervals).}
}
\examples{
xx <- c(1,2,1,2,1,2) #input factors
nlev <- rep(4,length(xx)) #number of levels for each factor
ntimes <- length(xx)/2 #base function is in 2D, so duplicate 3 times
camel6(xx,ntimes,nlev)
}
