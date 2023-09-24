unwgt_fRegress <- function(y, xfdlist, betalist,
                           y2cMap=NULL, SigmaE=NULL, returnMatrix=FALSE,
                           method=c('fRegress', 'model'),
                           sep='.', ...) {
  
  if (is.fdPar(y)) y <- y$fd
  
  arglist <- my_fRegressArgCheck(y, xfdlist, betalist)
  
  yfdobj   <- arglist$yfd  # the older version used yfdPar as the name.
  xfdlist  <- arglist$xfdlist
  betalist <- arglist$betalist
  
  p <- length(xfdlist)
  
  #  --------------------------------------------------------------------------
  #  branch depending on whether the dependent variable is functional or scalar
  #  --------------------------------------------------------------------------
  
  #  ----------------------------------------------------------------
  #           YFDOBJ is a functional data object
  #  ----------------------------------------------------------------
  
  #  extract dependent variable information
  ycoef     <- yfdobj$coefs
  ycoefdim  <- dim(ycoef)
  N         <- ycoefdim[2]
  ybasisobj <- yfdobj$basis
  rangeval  <- ybasisobj$rangeval
  ynbasis   <- ybasisobj$nbasis
  onesbasis <- create.constant.basis(rangeval)
  onesfd    <- fd(1,onesbasis)
  
  if (length(ycoefdim) > 2) stop("YFDOBJ from YFD is not univariate.")
  
  #  --------  set up the linear equations for the solution  -----------
  
  #  compute the total number of coefficients to be estimated
  
  ncoef <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      ncoefj     <- betafdParj$fd$basis$nbasis
      ncoef      <- ncoef + ncoefj
    }
  }
  
  Cmat <- matrix(0,ncoef,ncoef)
  Dmat <- rep(0,ncoef)
  
  #  ------------------------------------------------------------------------
  #  Compute the symmetric positive definite matrix CMAT and
  #  the column matrix DMAT.  CMAT contains weighted inner products of
  #  bases for each pair of terms plus, for lambda > 0, a roughness penalty
  #  matrix to ensure that the estimated coefficient functions will be smooth
  #  The weight vector is the point-wise product of the associated functional
  #  covariates.  
  #  Dmat contains for each covariate the weighted integral of the basis 
  #  functions, with the weight function being the covariate function
  #  pointwise multiplied the dependent variate yobj.
  #  The estimated coefficients functions are defined by the solution
  #  CMAT %*% COEF = DMAT.
  #  ------------------------------------------------------------------------
  
  #  loop through rows of CMAT
  
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      #  get jth beta basis
      betafdj    <- betafdParj$fd
      betabasisj <- betafdj$basis
      ncoefj     <- betabasisj$nbasis
      #  row indices of CMAT and DMAT to fill
      mj1    <- mj2 + 1
      mj2    <- mj2 + ncoefj
      indexj <- mj1:mj2
      #  compute right side of equation DMAT
      #  compute weight function for DMAT
      xfdj <- xfdlist[[j]]
      xyfdj <- xfdj*yfdobj
      wtfdj <- sum(xyfdj)
      #  Compute jth component of DMAT
      Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
      Dmat[indexj] <- Dmatj
      #  loop through columns of CMAT
      mk2 <- 0
      for (k in 1:j) {
        betafdPark <- betalist[[k]]
        if (betafdPark$estimate) {
          #  get the kth basis
          betafdk    <- betafdPark$fd
          betabasisk <- betafdk$basis
          ncoefk     <- betabasisk$nbasis
          #  column indices of CMAT to fill
          mk1 <- mk2 + 1
          mk2 <- mk2 + ncoefk
          indexk <- mk1:mk2
          #  set up weight function for CMAT component
          xfdk <- xfdlist[[k]]
          xxfdjk <- xfdj*xfdk
          wtfdjk <- sum(xxfdjk)
          #  compute the inner product
          Cmatjk <- inprod(betabasisj, betabasisk, 0, 0,
                           rangeval, wtfdjk)
          Cmat[indexj,indexk] <- Cmatjk
          Cmat[indexk,indexj] <- t(Cmatjk)
        }
      }
      #  attach penalty term to diagonal block if required
      lambdaj <- betafdParj$lambda
      if (lambdaj > 0) {
        Rmatj <- betafdParj$penmat
        if (is.null(Rmatj)) {
          Lfdj  <- betafdParj$Lfd
          Rmatj <- eval.penalty(betabasisj, Lfdj)
        }
        Cmat[indexj,indexj] <- Cmat[indexj,indexj] +
          lambdaj*Rmatj
      }
    }
  }
  
  #  ensure symmetry
  
  Cmat <- (Cmat+t(Cmat))/2
  
  Cmatinv <- solve(Cmat)
  
  betacoef <- Cmatinv %*% Dmat
  
  #  set up fdPar objects for reg. fns. in BETAESTLIST
  
  betaestlist <- betalist
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      betafdj <- betafdParj$fd
      ncoefj  <- betafdj$basis$nbasis
      mj1     <- mj2 + 1
      mj2     <- mj2 + ncoefj
      indexj  <- mj1:mj2
      coefj   <- betacoef[indexj]
      betafdj$coefs <- as.matrix(coefj)
      betafdParj$fd <- betafdj
    }
    betaestlist[[j]] <- betafdParj
  }
  
  #  set up fd objects for predicted values in YHATFDOBJ
  
  nfine   <- max(501,10*ynbasis+1)
  tfine   <- seq(rangeval[1], rangeval[2], len=nfine)
  yhatmat <- matrix(0,nfine,N)
  for (j in 1:p) {
    xfdj       <- xfdlist[[j]]
    xmatj      <- eval.fd(tfine, xfdj, 0, returnMatrix)
    betafdParj <- betaestlist[[j]]
    betafdj    <- betafdParj$fd
    betavecj   <- eval.fd(tfine, betafdj, 0, returnMatrix)
    yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
  }
  yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
  
  df <- NA
  
  #  -----------------------------------------------------------------------
  #        Compute pointwise standard errors of regression coefficients
  #               if both y2cMap and SigmaE are supplied.
  #        y2cMap is supplied by the smoothing of the data that defined
  #        the dependent variable.
  #        SigmaE has to be computed from a previous analysis of the data.
  #  -----------------------------------------------------------------------
  
  if (!(is.null(y2cMap) || is.null(SigmaE))) {
    
    #  check dimensions of y2cMap and SigmaE
    
    y2cdim = dim(y2cMap)
    if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) {
      stop("Dimensions of Y2CMAP not correct.")
    }
    
    ybasismat = eval.basis(tfine, ybasisobj, 0, returnMatrix)
    
    deltat    = tfine[2] - tfine[1]
    
    #  compute BASISPRODMAT
    
    basisprodmat = matrix(0,ncoef,ynbasis*N)
    
    mj2 = 0
    for (j in 1:p) {
      betafdParj = betalist[[j]]
      betabasisj = betafdParj$fd$basis
      ncoefj     = betabasisj$nbasis
      bbasismatj = eval.basis(tfine, betabasisj, 0, returnMatrix)
      xfdj       = xfdlist[[j]]
      tempj      = eval.fd(tfine, xfdj, 0, returnMatrix)
      #  row indices of BASISPRODMAT to fill
      mj1    = mj2 + 1
      mj2    = mj2 + ncoefj
      indexj = mj1:mj2
      #  inner products of beta basis and response basis
      #    weighted by covariate basis functions
      mk2 = 0
      for (k in 1:ynbasis) {
        #  row indices of BASISPRODMAT to fill
        mk1    = mk2 + 1
        mk2    = mk2 + N
        indexk = mk1:mk2
        tempk  = bbasismatj*ybasismat[,k]
        basisprodmat[indexj,indexk] =
          deltat*crossprod(tempk,tempj)
      }
    }
    
    #  compute variances of regression coefficient function values
    
    c2bMap    = solve(Cmat,basisprodmat)
    VarCoef   = y2cMap %*% SigmaE %*% t(y2cMap)
    CVariance = kronecker(VarCoef,diag(rep(1,N)))
    bvar      = c2bMap %*% CVariance %*% t(c2bMap)
    betastderrlist = vector("list", p)
    mj2 = 0
    for (j in 1:p) {
      betafdParj = betalist[[j]]
      betabasisj = betafdParj$fd$basis
      ncoefj     = betabasisj$nbasis
      mj1 	     = mj2 + 1
      mj2 	     = mj2 + ncoefj
      indexj     = mj1:mj2
      bbasismat  = eval.basis(tfine, betabasisj, 0, returnMatrix)
      bvarj      = bvar[indexj,indexj]
      bstderrj   = sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj = smooth.basis(tfine, bstderrj, betabasisj)$fd
      betastderrlist[[j]] = bstderrfdj
    }
  } else {
    betastderrlist = NULL
    bvar           = NULL
    c2bMap         = NULL
  }
  
  #  -------------------------------------------------------------------
  #                       Set up output list object
  #  -------------------------------------------------------------------
  
  fRegressList <-
    list(yfdobj         = yfdobj,
         xfdlist        = xfdlist,
         betalist       = betalist,
         betaestlist    = betaestlist,
         yhatfdobj      = yhatfdobj,
         Cmat           = Cmat,
         Dmat           = Dmat,
         Cmatinv        = Cmatinv,
         df             = df,
         y2cMap         = y2cMap,
         SigmaE         = SigmaE,
         betastderrlist = betastderrlist,
         bvar           = bvar,
         c2bMap         = c2bMap)
  
  
  return(fRegressList)
  
}

## Auxiliary functions ---------------------------------------------------------

my_fRegressArgCheck <- function(yfd, xfdlist, betalist) 
{
  #  FREGRESS_ARGCHECK checks the first three arguments for the functions
  #  for function regression, including FREGRESS.
  
  #  --------------------  Check classes of arguments  --------------------
  
  #  check that YFD is of class either 'fd' or 'numeric' and compute sample size N
  
  if (!(is.fdPar(yfd) || is.fd(yfd) || is.numeric(yfd) || is.matrix(yfd))) stop(
    "First argument is not of class 'fdPar', 'fd', 'numeric' or 'matrix'.")
  
  if (is.fdPar(yfd)) yfd <- yfd$fd
  
  if (inherits(yfd, "fd")) {
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
  } else {
    N <- length(yfd)
  } 
  
  #  check that xfdlist is a list object and compute number of covariates p
  
  #  check XFDLIST
  
  if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)
  
  if (!inherits(xfdlist, "list")) stop(
    "Argument XFDLIST is not a list object.")
  
  #  get number of independent variables p
  
  p <- length(xfdlist)
  
  #  check BETALIST
  
  if (inherits(betalist, "fd")) betalist <- list(betalist)
  
  if (!inherits(betalist, "list")) stop(
    "Argument BETALIST is not a list object.")
  
  if (length(betalist) != p)  {
    cat(paste("\nNumber of regression coefficients does not match\n",
              "number of independent variables."))
    stop("")
  }
  
  #  extract the range if YFD is functional
  
  if (inherits(yfd, "fd")) {
    rangeval <- yfd$basis$rangeval
  } else {
    rangeval = c(0,1)
  }
  
  #  --------------------  check contents of XFDLIST  -------------------
  
  #  If the object is a vector of length N,
  #  it is converted to a functional data object with a
  #  constant basis
  
  onebasis <- create.constant.basis(rangeval)
  onesfd   <- fd(1,onebasis)
  
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xcoef <- xfdj$coefs
      if (length(dim(xcoef)) > 2) stop(
        paste("Covariate",j,"is not univariate."))
      #  check size of coefficient array
      Nj <- dim(xcoef)[2]
      if (Nj != N) {
        print(
          paste("Incorrect number of replications in XFDLIST",
                "for covariate",j))
        xerror = TRUE
      }
    } 
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != N && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
    } 
    if (!(inherits(xfdlist[[j]], "fd"     ) || 
          inherits(xfdlist[[j]], "numeric") ||
          inherits(xfdlist[[j]], "matrix" ))) {
      print(paste("XFDLIST[[",j,"]] is not an FD or numeric or matrix object."))
      xerror = TRUE
    }
  }
  
  #  --------------------  check contents of BETALIST  -------------------
  
  berror <- FALSE
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (inherits(betafdParj, "fd") || inherits(betafdParj, "basisfd")) {
      betafdParj    <- fdPar(betafdParj)
      betalist[[j]] <- betafdParj
    }
    if (!inherits(betafdParj, "fdPar")) {
      print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
      berror <- TRUE
    }
  }
  
  if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")
  
  #  ---------------------  return the argument list  --------------------
  
  # The older versions of fda package used yfdPar as the name for the first member.
  
  return(list(yfd=yfd, xfdlist=xfdlist, betalist=betalist))
  
}
