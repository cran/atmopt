atm.init <- function(nfact, nlev){

  #initialize ATM object
  des.atm <- matrix(NA,nrow=0,ncol=nfact)
  des.all <- matrix(NA,nrow=0,ncol=nfact)
  obs.atm <- c()
  obs.all <- c()
  cur.lev <- vector("list",nfact) #list of current levels
  for (i in 1:nfact){
    cur.lev[[i]] <- 1:(nlev[i])
  }

  #return
  fit <- vector("list")
  fit$des.atm <- des.atm
  fit$des.all <- des.all
  fit$obs.atm <- obs.atm
  fit$obs.all <- obs.all
  fit$nfact <- nfact
  fit$nlev <- nlev
  fit$cur.lev <- cur.lev
  fit$nrem <- 0 #levels removed

  return(fit)

}

atm.nextpts <- function(atm.obj,reps=NULL){

  nfactor <- atm.obj$nfact
  nlev <- atm.obj$nlev
  levs.cur <- atm.obj$cur.lev
  if (is.null(reps)){
    reps <- 1
  }

  #Sample new OA and randomize
  OA <- oa.design(nfactors=nfactor,nlevels=nlev-atm.obj$nrem,replications=reps)
  des.mat <- matrix(NA,nrow=length(OA$A),ncol=nfactor)
  for (i in 1:nrow(des.mat)){
    for (j in 1:ncol(des.mat)){
      des.mat[i,j] <- OA[[j]][i]
    }
  }
  rnd.des <- rand.perm(des.mat)
  # print("Got here")

  #Return points
  new.pts <- trans.lvl(rnd.des,levs.cur)
  # print("Got here")

}

atm.addpts <- function(atm.obj,des.new,obs.new){

  fit <- atm.obj

  #Add new design points and observations
  fit$des.atm <- rbind(fit$des.atm,des.new)
  fit$des.all <- rbind(fit$des.all,des.new)
  fit$obs.atm <- c(fit$obs.atm,obs.new)
  fit$obs.all <- c(fit$obs.all,obs.new)

  return(fit)
}

atm.predict <- function(atm.obj,alphas=NULL,ntimes=1,nsub=100,
                        prob.am=0.5,prob.pw=1.0,reps=NULL){

  if (is.null(reps)){
    reps <- 1
  }

  #Tune alphas
  if (is.null(alphas)){
    alphas <- tune.alpha(atm.obj$des.atm,c(atm.obj$obs.atm),nlevels=atm.obj$nlev-atm.obj$nrem,rp=reps,
                         des.all=atm.obj$des.all,nboot=ntimes,subsamp=nsub,
                         prob.am=prob.am,prob.pw=prob.pw)
  }

  #Choose index for prediction
  fn.lst.atm <- cte.fnl(alphas)
  mm <- marg.min(atm.obj$des.atm,atm.obj$obs.atm,fn.lst.atm)
  if (all(alphas==0)){
    #if alpha = 0, then take the best observed solution
    idx.opt <- atm.obj$des.all[which.min(atm.obj$obs.all),]
  }else{#otherwise take ATM prediction
    idx.opt <- mm$margmin
  }

  #Update ATM object
  atm.obj$mm <- mm
  atm.obj$idx.opt <- idx.opt

  return(atm.obj)
}

atm.remlev <- function(atm.obj){

  tmp <- remove.levels(atm.obj$des.atm,atm.obj$obs.atm,atm.obj$cur.lev,(atm.obj$mm)$margmax)
  atm.obj$cur.lev <- tmp$lst
  atm.obj$des.atm <- tmp$des
  atm.obj$obs.atm <- tmp$obs
  atm.obj$nrem <- atm.obj$nrem + 1

  return(atm.obj)

}
