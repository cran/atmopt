marg.min <- function(des,obs,fn.lst,nmax=NULL){
  # Computes an estimate of the minimum using marginal information
  stg <- rep(NA,ncol(des)) # for each factor
  if (is.null(nmax)){
    stm <- rep(NA,ncol(des))
  }else{
    stm <- vector("list",ncol(des))
  }
  smy <- vector("list",ncol(des))
  levs <- vector("list",ncol(des))

  #Get unique levels for each factor
  for (i in 1:ncol(des)){
    levs[[i]] <- unique(des[,i])
    smy[[i]] <- rep(NA,length(levs[[i]]))
  }

  #Find the best marginal level
  for (i in 1:ncol(des)){
    # smy <- rep(NA,max(des[,i]))
    fn <- fn.lst[[i]]$fn
    for (j in 1:length(levs[[i]])){# for each level of factor i
      smy[[i]][j] <- fn(obs[which(des[,i]==levs[[i]][j])])
    }
    stg[i] <- levs[[i]][which.min(smy[[i]])]
    if (is.null(nmax)){
      stm[i] <- levs[[i]][which.max(smy[[i]])]
    }else{
      ord <- order(smy[[i]],decreasing=TRUE)
      stm[[i]] <- levs[[i]][ord[1:nmax[i]]] #take largest indices to remove
    }
  }

  #if it's minimum, just do which.min (marginal mins fail if same obs)
  ind.almn <- TRUE;
  for (i in 1:ncol(des)){
    ind.almn <- ind.almn & (fn.lst[[i]]$alpha==0)
  }
  if (ind.almn){#if all zeros
    # print("Defaulting to PW...")
    margmin = des[which.min(obs),]
  }

  return(list(margmin=stg,margmax=stm,margsum=smy))
}

remove.levels <- function(des,obs,lst.lev,mx){
  lst.ret <- lst.lev
  for (i in 1:length(lst.lev)){
    vec <- lst.lev[[i]]
    if (length(vec)>1){
      if (!is.list(mx)){
        lst.ret[[i]] <- vec[which(vec!=mx[i])]
      }else{
        lst.ret[[i]] <- setdiff(vec,match(mx[[i]],vec))
      }
    }
  }

  # des.ret <- des
  idx.rm <- c() #indices to remove
  for (i in 1:ncol(des)){
    vec <- lst.lev[[i]]
    if (length(vec)>1){
      if (!is.list(mx)){
        idx.rm <- c(idx.rm, which(des[,i]==mx[i]))
      }else{
        idx.tmp <- match(mx[[i]],des[,i])
        idx.rm <- c(idx.rm, idx.tmp[!is.na(idx.tmp)])
      }
    }
  }
  idx.rm <- unique(idx.rm) #take unique indices
  des.ret <- des[-idx.rm,]
  obs.ret <- obs[-idx.rm]

  return(list(des=des.ret,obs=obs.ret,lst=lst.ret))
}
