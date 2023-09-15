## Weighted penalized smoothing on cubic B-spline basis
wt_bsplinesmoothing <- function(curves, wgts.obs, t.points, breaks, lambda=NULL, set.cb = FALSE)
{
  ## Utilities
  n      <- dim(curves)[2]
  
  ## Define cubic B-spline basis
  basis  <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
  
  L      <- basis$nbasis
  wgts   <- wgts.obs
  
  if(is.null(lambda))
  {
    esp        <- seq(0,7, by=1)
    lambda.vec <- 10^-esp
    gcv.mat    <- matrix(data=0, nrow=dim(curves)[2], ncol=length(lambda.vec))
    for(j in 1:length(lambda.vec))
    {
      functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
      for(i in 1:n)
      {
        curve            <- curves[,i]
        wtvec            <- wgts[,i]
        curve.s          <- smooth.basis(t.points, curve, functionalPar, wtvec=wtvec)
        gcv.mat[i,j]     <- curve.s$gcv
      }
    }
    gcv <- colSums(gcv.mat)
    lambda.opt <- lambda.vec[which(gcv == min(gcv))]
    
  } else {
    lambda.opt  <- lambda
  }
  
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
  
  # Smooth curves
  curves.scoef  <- matrix(data=0, nrow=L, ncol=n)
  N <- length(t.points)
  y2cmaps <- array(data=0, dim=c(L,N,n))
  for(i in 1:n)
  {
    curve            <- curves[,i]
    wtvec            <- wgts[,i]
    curve.s          <- smooth.basis(t.points, curve, fPar, wtvec=wtvec)
    curves.scoef[,i] <- curve.s$fd$coefs
    
    if(set.cb)
    {
      y2cmaps[,,i] <- curve.s$y2cMap
    }
  }
  
  curves.s          <- smooth.basis(t.points, curves, fPar)
  curves.s$fd$coefs <- curves.scoef
  
  ## Output
  out        <- list(curves.s$fd, wgts.obs, fPar, lambda.opt, y2cmaps)
  names(out) <- list('curves.fd', 'wgts.obs', 'fPar', 'lambda.opt', 'y2cmaps')
  
  return(out)
  
}
