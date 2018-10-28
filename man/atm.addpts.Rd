\name{atm.addpts}
\alias{atm.addpts}
\docType{package}
\title{
  Add new data for ATM
}
\description{
  \code{atm.addpts} adds new data into an ATM object.
}
\usage{
  atm.addpts(atm.obj,des.new,obs.new)
}
\arguments{
  \item{atm.obj}{Current ATM object.}
  \item{des.new}{Design matrix for new evaluations.}
  \item{obs.new}{Observations for new evaluations.}
}
\examples{
  \dontrun{
####################################################
# Example 1: detpep10exp (9-D)
####################################################
nfact <- 9 #number of factors
ntimes <- floor(nfact/3) #number of "repeats" for detpep10exp
lev <- 4 #number of levels
nlev <- rep(lev,nfact) #number of levels for each factor
nelim <- 3 #number of level eliminations
fn <- function(xx){detpep10exp(xx,ntimes,nlev)} #objective to minimize (assumed expensive)

#initialize objects
# (predicts & removes levels based on tuned ATM percentages)
fit.atm <- atm.init(nfact,nlev)
#initialize sel.min object
# (predicts minimum using smallest observed value & removes levels with largest minima)
fit.min <- atm.init(nfact,nlev)

#Run for nelim eliminations:
res.atm <- rep(NA,nelim) #for ATM results
res.min <- rep(NA,nelim) #for sel.min results
for (i in 1:nelim){

  # ATM updates:
  new.des <- atm.nextpts(fit.atm) #get design points
  new.obs <- apply(new.des,1,fn) #sample function
  fit.atm <- atm.addpts(fit.atm,new.des,new.obs) #add data to object
  fit.atm <- atm.predict(fit.atm) #predict minimum setting
  idx.atm <- fit.atm$idx.opt
  res.atm[i] <- fn(idx.atm)
  fit.atm <- atm.remlev(fit.atm) #removes worst performing level

  # sel.min updates:
  new.des <- atm.nextpts(fit.min) #get design points
  new.obs <- apply(new.des,1,fn) #sample function
  fit.min <- atm.addpts(fit.min,new.des,new.obs) #add data to object
  fit.min <- atm.predict(fit.min, alphas=rep(0,nfact)) #find setting with smallest observation
  idx.min <- fit.min$idx.opt
  res.min[i] <- fn(idx.min)
  #check: min(fit.min$obs.all)
  fit.min <- atm.remlev(fit.min) #removes worst performing level

}

res.atm
res.min

#conclusion: ATM finds better solutions by learning & exploiting additive structure

####################################################
# Example 2: camel6 (24-D)
####################################################
nfact <- 24 #number of factors
ntimes <- floor(nfact/2.0) #number of "repeats" for camel6
lev <- 4
nlev <- rep(lev,nfact) #number of levels for each factor
nelim <- 3 #number of level eliminations
fn <- function(xx){camel6(xx,ntimes,nlev)} #objective to minimize (assumed expensive)

#initialize objects
# (predicts & removes levels based on tuned ATM percentages)
fit.atm <- atm.init(nfact,nlev)
#initialize sel.min object
# (predicts minimum using smallest observed value & removes levels with largest minima)
fit.min <- atm.init(nfact,nlev)

#Run for nelim eliminations:
res.atm <- rep(NA,nelim) #for ATM results
res.min <- rep(NA,nelim) #for sel.min results
for (i in 1:nelim){

  # ATM updates:
  new.des <- atm.nextpts(fit.atm) #get design points
  new.obs <- apply(new.des,1,fn) #sample function
  fit.atm <- atm.addpts(fit.atm,new.des,new.obs) #add data to object
  fit.atm <- atm.predict(fit.atm) #predict minimum setting
  idx.atm <- fit.atm$idx.opt
  res.atm[i] <- fn(idx.atm)
  fit.atm <- atm.remlev(fit.atm) #removes worst performing level

  # sel.min updates:
  new.des <- atm.nextpts(fit.min) #get design points
  new.obs <- apply(new.des,1,fn) #sample function
  fit.min <- atm.addpts(fit.min,new.des,new.obs) #add data to object
  fit.min <- atm.predict(fit.min, alphas=rep(0,nfact)) #find setting with smallest observation
  idx.min <- fit.min$idx.opt
  res.min[i] <- fn(idx.min)
  #check: min(fit.min$obs.all)
  fit.min <- atm.remlev(fit.min) #removes worst performing level

}

res.atm
res.min

#conclusion: ATM finds better solutions by learning & exploiting additive structure
}
}
