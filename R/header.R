rand.perm <- function(des){
  # Randomly permutes the labeling for a design
  newdes <- des
  for (i in 1:ncol(des)){
    mx <- max(des[,i])
    prm <- sample.int(mx)
    for (j in 1:nrow(des)){
      newdes[j,i] <- prm[des[j,i]]
    }
  }
  return(newdes)
}

trans.uh <- function(des,nlevels,levs=NULL){
  #des - design in original levels
  if (length(nlevels)==1){
    nlevels <- rep(nlevels,ncol(des))
  }
  des.mat <- des
  for (i in 1:ncol(des)){
    for ( j in 1:nrow(des)){
      if (is.null(levs)){
        #if levels not specified, take midpoint
        des.mat[j,i] <- (des[j,i]-0.5)/nlevels[i]
      }else{
        #else take given levels
        des.mat[j,i] <- levs[i,des[j,i]]
      }
    }
  }
  return(des.mat)
}

trans.lvl <- function(des,lst.lev){
  des.ret <- des
  for (i in 1:nrow(des)){
    for (j in 1:ncol(des)){
      des.ret[i,j] <- lst.lev[[j]][des[i,j]]
    }
  }
  return(des.ret)
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

# hierNet.path.mod <- function (x, y, lamlist = NULL, delta = 1e-08, minlam = NULL,
#           maxlam = NULL, nlam = 20, flmax = 1, flmin = 0.01, diagonal = TRUE,
#           strong = FALSE, aa = NULL, zz = NULL, stand.main = TRUE,
#           stand.int = FALSE, rho = nrow(x), niter = 100, sym.eps = 0.001,
#           step = 1, maxiter = 2000, backtrack = 0.2, tol = 1e-05, trace = 0)
# {
#   this.call <- match.call()
#   x <- scale(x, center = TRUE, scale = stand.main)
#   mx <- attr(x, "scaled:center")
#   sx <- attr(x, "scaled:scale")
#   my <- mean(y)
#   y <- y - my
#   if (is.null(maxlam)) {
#     if (!is.null(minlam))
#       stop("Cannot have maxlam=NULL if minlam is non-null.")
#     maxlam <- max(abs(t(x) %*% y)) * flmax
#     minlam <- maxlam * flmin
#   }
#   if (is.null(minlam))
#     minlam <- maxlam * flmin
#   if (is.null(lamlist))
#     lamlist <- exp(seq(log(maxlam), log(minlam), length = nlam))
#   nlam <- length(lamlist)
#   if (is.null(zz))
#     zz <- compute.interactions.c(x, diagonal = diagonal)
#   else stopifnot(is.matrix(zz))
#   zz <- scale(zz, center = TRUE, scale = stand.int)
#   mzz <- attr(zz, "scaled:center")
#   szz <- attr(zz, "scaled:scale")
#   zz <- as.numeric(zz)
#   p <- ncol(x)
#   cp2 <- choose(p, 2)
#   bp <- bn <- matrix(NA, nrow = p, ncol = nlam)
#   th <- array(NA, c(p, p, nlam))
#   obj <- rep(NA, nlam)
#   aa <- NULL
#   for (i in seq(nlam)) {
#     cat(c("i,lam=", i, round(lamlist[i], 2)), fill = TRUE)
#     aa <- hierNet(x, y, lam = lamlist[i], delta = delta,
#                   strong = strong, diagonal = diagonal, aa = aa, zz = zz,
#                   stand.main = FALSE, stand.int = FALSE, rho = rho,
#                   niter = niter, sym.eps = sym.eps, step = step, maxiter = maxiter,
#                   backtrack = backtrack, tol = tol, trace = trace)
#     bp[, i] <- aa$bp
#     bn[, i] <- aa$bn
#     th[, , i] <- aa$th
#     obj[i] <- aa$obj
#   }
#   dimnames(bp) <- dimnames(bn) <- list(as.character(1:p), NULL)
#   dimnames(th) <- list(as.character(1:p), as.character(1:p),
#                        NULL)
#   out <- list(bp = bp, bn = bn, th = th, obj = obj, lamlist = lamlist,
#               delta = delta, mx = mx, sx = sx, mzz = mzz, szz = szz,
#               my = my, type = "gaussian", diagonal = diagonal, strong = strong,
#               step = step, maxiter = maxiter, backtrack = backtrack,
#               tol = tol, call = this.call)
#   if (strong) {
#     out$rho <- rho
#     out$niter <- niter
#     out$sym.eps <- sym.eps
#   }
#   class(out) <- "hierNet.path"
#   out
# }
