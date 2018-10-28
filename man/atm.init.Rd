\name{atm.init}
\alias{atm.init}
\docType{package}
\title{
  Initializing ATM object
}
\description{
  \code{atm.init} initialize the ATM object to use for optimization.
}
\usage{
  atm.init(nfact, nlev)
}
\arguments{
  \item{nfact}{Number of factors to optimize.}
  \item{nlev}{A vector containing the number of levels for each factor.}
}
\examples{
nfact <- 9 #number of factors
lev <- 4
nlev <- rep(lev,nfact) #number of levels for each factor
fit <- atm.init(nfact,nlev) #initialize ATM object
}
