tune.alpha <- function(des,obs,nlevels,rp,des.all=des,
                       nboot=25,subsamp=100,prob.pw=0.5,prob.am=0.5){
  nfactor <- ncol(des)

  #Fit interactions model
  fac.mat <- data.frame(des)
  fac.mat <- lapply(fac.mat,factor)
  mod.mtx <- model.matrix(~.,fac.mat)
  mod.mtx <- mod.mtx[,-1] #remove intercept
  fit <- hierNet.path(mod.mtx,obs,diagonal=FALSE,stand.main=FALSE,flmin=1e-4)
  # fit <- hierNet.path.mod(mod.mtx,obs,diagonal=FALSE,stand.main=FALSE,flmax=0.01,flmin=1e-4)
  # fit <- glmnet(mod.mtx,obs)
  nf.flg <- TRUE
  nfcur <- 10
  while (nf.flg){
    tryCatch(
      {
        fitcv <- hierNet.cv(fit,des,obs,nfolds=nfcur);
        # fitcv <- cv.glmnet(des,obs,nfold=nfcur,lambda=fit$lambda)
        nf.flg <- FALSE
      },
      error=function(cond){},
      finally={nfcur <- ceiling(nfcur / 2.0)}
    )
  }
  ind <- which.min(abs(fit$lamlist-fitcv$lamhat)) #set best lambda

  #Compute optimal settings for each alpha
  if (all(nlevels==nlevels[1])){
    alpha.vec <- seq(from=0,to=1,length.out=nlevels[1]+1)
    if (is.na(subsamp)){
      alpha.mat <- permutations(length(alpha.vec),nfactor,v=alpha.vec,repeats.allowed=T)
    }else{
      # alpha.mat <- permutations(length(alpha.vec),nfactor,v=alpha.vec,repeats.allowed=T)
      if (prod(nlevels)>subsamp){
        alpha.mat <- matrix(sample(alpha.vec,nfactor*subsamp,replace=T),ncol=nfactor)
      }else{
        alpha.mat <- permutations(length(alpha.vec),nfactor,v=alpha.vec,repeats.allowed=T)
      }
    }
  }else{#uneven levels
    alpha.mat <- matrix(NA,nrow=subsamp,ncol=nfactor)
    for (i in 1:nfactor){
      alpha.vec <- seq(from=0,to=1,length.out=nlevels[i]+1)
      alpha.mat[,i] <- sample(alpha.vec,subsamp,replace=T)
    }
  }
  alpha.mat <- rbind(rep(0.0,nfactor),alpha.mat) #Add PW
  alpha.mat <- rbind(rep(1.0,nfactor),alpha.mat) #Add AM

  #Bootstrap from fitted model
  idx.mat <- matrix(NA,nrow=nboot,ncol=nrow(alpha.mat))
  for (mm in 1:nboot){
    # generate random OA for testing
    OA <- oa.design(nfactors=nfactor,nlevels=nlevels,replications=rp)
    des.rnd <- matrix(NA,nrow=length(OA$A),ncol=nfactor)
    for (i in 1:nrow(des.rnd)){
      for (j in 1:ncol(des.rnd)){
        des.rnd[i,j] <- OA[[j]][i]
      }
    }
    des.rnd <- rand.perm(des.rnd) #random level permutation
    fac.mat <- data.frame(des.rnd) #get indicator matrix
    fac.mat <- lapply(fac.mat,factor)
    mod.mtx <- model.matrix(~.,fac.mat)
    mod.mtx <- mod.mtx[,-1] #remove intercept
    obs.rnd <- predict(fit,mod.mtx,lam=fitcv$lamhat) #predict
    obs.rnd <- obs.rnd[,ind]

    # compute predictor
    ll.vec <- matrix(NA,nrow=nrow(alpha.mat),ncol=nfactor)
    for (i in 1:nrow(alpha.mat)){
      alphas <- alpha.mat[i,]
      fn.lst <- cte.fnl(alphas)
      ll.vec[i,] <- marg.min(des.rnd,obs.rnd,fn.lst)$margmin
      # ll.vec[i,] <- marg.min(des,obs,fn.lst)$margmin
    }
    new.mat <- data.frame(ll.vec) #compute model matrix for dummy variables
    new.lst <- vector("list",ncol(new.mat))
    lev.strs <- lapply(fac.mat,levels)
    for (j in 1:ncol(new.mat)){
      new.mat[,j] <- factor(new.mat[,j],levels=lev.strs[[j]])
    }
    new.mtx <- model.matrix(~.,new.mat,contrasts.arg=lapply(new.mat,contrasts))
    new.mtx <- new.mtx[,-1] #remove intercept

    # #Refit least squares
    # actme <- (fit$bp[,ind]!=0)+(fit$bn[,ind]!=0)
    # actint <- which(fit$th[,,ind]!=0,arr.ind=T) #refit with least squares
    # newxx <- cbind(mod.mtx[,actme], apply(actint,1,function(xx){mod.mtx[,xx[1]]*mod.mtx[,xx[2]]}))

    #Tune alphas
    pp <- predict(fit,new.mtx,lam=fitcv$lamhat)
    pp <- pp[,ind]
    idx <- which(pp==min(pp))
    idx.mat[mm,] <- (pp==min(pp))
  }
  idx <- which(colSums(idx.mat)==max(colSums(idx.mat)))

  idx.l <- idx #indices left
  #Remove indices which we already observed, unless it's the smallest
  for (ii in 1:length(idx)){
    fn.lst.chk <- cte.fnl(alpha.mat[idx[ii],])
    mm.chk <- marg.min(des,obs,fn.lst.chk)
    idx.chk <- mm.chk$margmin
    if (any(apply(des.all,1,function(xx){all(xx==idx.chk)}))){
      #remove this as a potential point to explore (this is incorporated later in PW)
      idx.l <- idx.l[-(idx.l==idx[ii])]
    }
  }
  idx <- idx.l

  #Randomize alphas if ties
  if (length(idx)==0){
    #if nothing left, that means everything was observed from data. then PW
    idx.sel <- 2
  }else if (length(idx)==1){
    #if only one choice, then take it!
    idx.sel <- idx
  }else if (2%in%idx){
    #If PW in choice, then randomly choose PW with prob.pw
    if (runif(1) <= prob.pw){
      idx.sel <- 2
    }else{
      idx.sel <- sample(idx,1)
    }
  }else if(1%in%idx){
    #If AM in choice, then randomly choose AM with prob.am
    if (runif(1) <= prob.am){
      idx.sel <- 1
    }else{
      idx.sel <- sample(idx,1)
    }
  }else{
    idx.sel <- sample(idx,1)
  }
  alphas <- alpha.mat[idx.sel,]

  return(alphas)
}
