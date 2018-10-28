\name{atm.predict}
\alias{atm.predict}
\docType{package}
\title{
  Predict the minimum setting for ATM
}
\description{
  \code{atm.init} predicts the minimum setting for an ATM object.
}
\usage{
  atm.predict(atm.obj,alphas=NULL,ntimes=1,nsub=100,prob.am=0.5,prob.pw=1.0,reps=NULL)
}
\arguments{
  \item{atm.obj}{Current ATM object.}
  \item{alphas}{A \code{p}-vector for ATM percentiles. \code{NULL} if tuned from data.}
  \item{ntimes}{Number of resamples for tuning ATM percentages.}
  \item{nsub}{Number of candidate percentiles to consider.}
  \item{prob.am}{In case of ties in percentage estimation, the probability of choosing marginal means (if optimal) for minimization.}
  \item{prob.pw}{In case of ties in percentage estimation, probability of picking-the-winner (if optimal) for minimization.}
  \item{reps}{Number of replications for internal OA in tuning ATM percentages.}
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
