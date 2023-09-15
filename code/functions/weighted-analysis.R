
## Workflow weighted analysis -----------------------------------------------------
workflow_weighted_analysis <- function(b,
                                       B,
                                       t.points,
                                       breaks,
                                       T_hp,
                                       curves,
                                       curves.true.fd,
                                       xlist,
                                       method = c('Kraus', 'KLNoAl' ,'KLAl', NULL),
                                       fix.par = 10,
                                       wgts.flag = TRUE,
                                       wgts.recon.flag=FALSE)
{
  
  #'# Load methods ---------------------------------------
  source('functions/weighted-smoothing.R')
  source('functions/weighted-regression.R')
  source('functions/reconstruction.R')
  
  ## Utilities for the identification of the b-th batch
  n <- dim(curves)[2]
  p <- length(xlist)
  
  #'# Separate training and test set ---------------------
  n.test <- floor(n/B)
  idxs <- ((b-1)*n.test + 1):min((b*n.test),n)
  
  test <- rep(FALSE, n)
  test[idxs] <- TRUE
  train <- (!test)
  
  curves.test    <- curves[,test]
  curves.train   <- curves[,train]
  
  # Create xlist.test and xlist.train
  xlist.test  <- list()
  xlist.train <- list()
  
  for(i in 1:length(xlist))
  {
    if (inherits(xlist[[i]], "fd"))
    {
      xlist.test[[i]]  <- xlist[[i]]
      xlist.test[[i]]$coefs  <- xlist[[i]]$coefs[,test]
      
      xlist.train[[i]] <- xlist[[i]]
      xlist.train[[i]]$coefs <- xlist[[i]]$coefs[,train]
      
    } else if (inherits(xlist[[i]], "numeric")) {
      
      xlist.test[[i]]  <- xlist[[i]][test]
      xlist.train[[i]] <- xlist[[i]][train]
      
    } else if (inherits(xlist[[i]], "matrix" )) {
      
      xlist.test[[i]]  <- xlist[[i]][test,1]
      xlist.train[[i]] <- xlist[[i]][train,1]
    }
  }
  
  reconst_fcts   <- find_obs_inc(Y = curves.train)
  
  T_hp.test      <- T_hp[test]
  T_hp.train     <- T_hp[train]
  
  ## Reconstruction ----------------------------------------------------------
  ## Kraus method
  if(any(method == 'Kraus'))
  {
    reconstruction <- reconstructKraus(X_mat        = curves.train,
                                       alpha        = NULL,
                                       reconst_fcts = reconst_fcts)
    
    curves.recon <- reconstruction[['X_reconst_mat']]
    wgts.recon   <- reconstruction[['W_reconst_mat']]
  }
  
  ## Kneip-Liebl Yes Alignment
  if(any(method == 'KLAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    reconstruction <- reconstructKneipLiebl(Ly             = Y_list,
                                            Lu             = U_list,
                                            method         = 'Error>0_AlignYES_CEscores',
                                            K              = NULL,
                                            reconst_fcts   = reconst_fcts,
                                            nRegGrid       = NULL,
                                            maxbins        = NULL,
                                            progrbar       = FALSE)
    
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                           nrow = length(t.points), ncol = length(reconst_fcts))
    
    wgts.recon <- matrix(unlist(reconstruction[['W_reconst_list']]),
                         nrow = length(t.points), ncol = length(reconst_fcts))
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves.train[,reconst_fcts[ii]]),ii] <- curves.train[!is.na(curves.train[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    
  }
  
  ## Kneip-Liebl No Alignment
  if(any(method == 'KLNoAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    reconstruction <- reconstructKneipLiebl(Ly           = Y_list,
                                            Lu           = U_list,
                                            method       = 'Error>0_AlignNO_CEscores',
                                            K            = NULL,
                                            reconst_fcts = reconst_fcts,
                                            nRegGrid     = NULL,
                                            maxbins      = NULL,
                                            progrbar     = FALSE)
    
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                           nrow = length(t.points), ncol = length(reconst_fcts))
    
    wgts.recon <- matrix(unlist(reconstruction[['W_reconst_list']]),
                         nrow = length(t.points), ncol = length(reconst_fcts))
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves.train[,reconst_fcts[ii]]),ii] <- curves.train[!is.na(curves.train[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    
  }
  
  curves.train.rec <- curves.train
  curves.train.rec[,reconst_fcts] <- curves.recon
  
  # curve.rec.Kraus <- curves.train.rec[,reconst_fcts[1]]
  # curve.rec.KLAL <- curves.train.rec[,reconst_fcts[1]]
  # curve.rec.KLPC <- curves.train.rec[,reconst_fcts[1]]
  
  # pdf(file = paste0("Results/images/po-curve-and-reconstruction.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(t.points, curves.train[, reconst_fcts[1]], #type='l', lwd=3,
  #      col='darkblue', pch=1,
  #      main="(a) Partially observed curve", ylim=c(-3.5,-1), xlab="t", ylab="y(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(t.points[t.points>1.75], curves.train.rec[(t.points>1.75), reconst_fcts[1]],
  #       col='darkblue', pch=16)
  # grid()
  # dev.off()
  
  # pdf(file = paste0("Results/images/po-curve-and-reconstruction-revision.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(t.points, curve.rec.KLAL, #type='l', lwd=3,
  #      col='darkblue', pch=1,
  #      main="(a) Partially observed curve", ylim=c(-3.5,-1), xlab="t", ylab="y(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(t.points[t.points>1.75], curve.rec.Kraus[(t.points>1.75)],
  #        col='orange', pch=17)
  # points(t.points[t.points>1.75], curve.rec.KLAL[(t.points>1.75)],
  #        col='darkblue', pch=19)
  # points(t.points[t.points>1.75], curve.rec.KLPC[(t.points>1.75)],
  #        col='forestgreen', pch=15)
  # legend(0,-2.5, legend=c('Kraus', 'KL-AL','KL-PC'), col=c('orange','darkblue', 'forestgreen'),
  #        pch=c(17,19,15), cex=1.6)
  # grid()
  # dev.off()
  
  #' Representation of one curve in its observed part and in its missing part
  #' (answer to reviewer)
  #' 
  # true.train <- curves.true[,train]
  # curve.true <- true.train[,reconst_fcts[1]]
  # 
  # pdf(file = paste0("Results/images/po-curve-and-missingpart.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(t.points, curve.true, #type='l', lwd=3,
  #      col='darkblue', pch=1,
  #      main="(b) True curve over the whole domain", ylim=c(-3.5,-1), xlab="t", ylab="y(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # grid()
  # dev.off()
  
  
  ## Definition of the weights -------------------------------------------------
  if(is.numeric(fix.par))
  {
    if(wgts.recon.flag){
      wgts.obs <- matrix(data=1, nrow=dim(curves.train)[1], ncol=dim(curves.train)[2])
      wgts.obs[,reconst_fcts] <- wgts.recon
      basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=2)
      wgts.fd  <- smooth.basis(t.points, wgts.obs, basis)$fd
      wgt <- list('wgts.obs'= wgts.obs,
                  'wgts.fd' = wgts.fd)
      
      # wgts.obs.Kraus <- wgts.obs[,reconst_fcts[1]]
      # wgts.obs.KLAL <- wgts.obs[,reconst_fcts[1]]
      # wgts.obs.KLPC <- wgts.obs[,reconst_fcts[1]]
      # save(wgts.obs.Kraus, wgts.obs.KLAL, wgts.obs.KLPC, file="Results/images/plot-recwgts-info.RData")
      # load("Results/images/plot-recwgts-info.RData")
      # pdf(file = paste0("Results/images/rec-weight-revision.pdf"), width = 8, height = 5)
      # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
      # plot(t.points, wgts.obs.KLAL, #type='l', lwd=3,
      #      col='darkblue', pch=19,
      #      main="(b) Reconstruction-driven weight", ylim=c(0,1.1), xlab="t", ylab="w(t)",
      #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
      # points(t.points, wgts.obs.Kraus, col='orange', pch=17)
      # legend(0,0.35, legend=c('Kraus', 'KL-AL/KL-PC'),
      #        col=c('orange','darkblue'),
      #        pch=c(17,19), cex=1.6)
      # grid()
      # dev.off()
      
    } else {
      wgt       <- create_weights(curves.rec    = curves.train.rec,
                                  t.points      = t.points,
                                  breaks        = breaks,
                                  fix.par       = fix.par,
                                  reconst_fcts  = reconst_fcts,
                                  Thp           = T_hp.train)
    }
    
    
    # plot(wgt$wgts.fd, ylim=c(0,1))
    # title(main=paste0('Parameter ', fix.par))
  } else {
    wgt       <- create_zero_weights(curves.rec    = curves.train.rec,
                                     t.points      = t.points,
                                     breaks        = breaks,
                                     reconst_fcts  = reconst_fcts,
                                     Thp           = T_hp.train)
  }
  
  ## Smoothing -----------------------------------------------------------------
  if(wgts.flag) # && is.numeric(fix.par)
  {
    smth      <- wt_bsplinesmoothing(curves   = curves.train.rec,
                                     wgts.obs = wgt$wgts.obs,
                                     t.points = t.points,
                                     breaks   = breaks,
                                     lambda   = NULL)
    curves.fd <- smth$curves.fd
    fPar  <- smth$fPar
  } else {
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    esp        <- seq(-7,1, by=1)
    lambda.vec <- sort(10^esp)
    gcv.vec    <- numeric(length(lambda.vec))
    for(j in 1:length(lambda.vec))
    {
      fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
      gcv.vec[j] <- sum(smooth.basis(t.points, curves.train.rec, fPar)$gcv)/n
    }
    lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
    curves.fd <- smooth.basis(t.points, curves.train.rec, fPar)$fd
  }
  
  ## B-list
  blist <- list(fPar,fPar,fPar)
  
  ## Regression and beta estimation --------------------------------------------
  if(wgts.flag)
  {
    mod <- weighted_fRegress(y            = curves.fd,
                             xfdlist      = xlist.train,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
  } else {
    mod <- fRegress(y            = curves.fd,
                    xfdlist      = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
  }
  
  beta_estimates <- mod$betaestlist
  ## Prepare the list of functional covariates which we use for prediction on
  ## test set
  onebasis <- create.constant.basis(range(t.points))
  onesfd   <- fd(1,onebasis)
  xfdlist.test  <- xlist.test
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist.test[[j]]
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != n.test && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist.test[[j]] <- fd(matrix(xfdj,1,n.test), onebasis)
    }
  }
  
  ## Error evaluation ----------------------------------------------------------
  curves.hat   <- predict_fRegress(mod          = mod,
                                   xlist        = xfdlist.test,
                                   t.points     = t.points)
  curves.hat.vals <- eval.fd(t.points, curves.hat)
  
  
  ## Evaluation of the MSE as integral mean of the squared functional errors
  curves.true.fd.test <- curves.true.fd
  curves.true.fd.test$coefs <- as.matrix(curves.true.fd$coefs[,test], ncol=n.test)
  
  MSE <- eval_MSE_functional(curves     = curves.true.fd.test,
                             curves.hat = curves.hat,
                             t.points   = t.points)
  
  
  ## Output --------------------------------------------------------------------
  out.list <- list(MSE                = MSE,
                   beta_estimates     = beta_estimates)
  
  return(out.list)
  
}


## Auxiliary functions ---------------------------------------------------------

find_obs_inc <- function(Y)
{
  n            <- dim(Y)[2]
  incomplete   <- c()
  for(i in 1:n)
  {
    incomplete[i] <- any(is.na(Y[,i]))
  }
  reconst_fcts <- which(incomplete==TRUE)
  return(reconst_fcts)
}


#### Definition of logistic weights
create_weights <- function(curves.rec, t.points, breaks, fix.par, reconst_fcts, Thp)
{
  ## Utilities
  n <- dim(curves.rec)[2]
  
  ## Creating the weights
  #tt <- seq(range(t.points)[1], range(t.points)[2], by=0.01)
  tt <- t.points
  
  T.last <- range(t.points)[2]
  wgts.obs <- matrix(data=1, nrow=length(tt), ncol=n)
  
  for(i in 1:length(reconst_fcts)) # i = 1
  {
    t.max    <- Thp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    T.max    <- t.points[imax]
    
    scale    <- fix.par*sd(curves.rec[imax,-reconst_fcts])
    loc      <- T.max + (T.last - T.max)/2
    
    # if(set.log)
    # {
    #   loc    <- 10^T.max + loc.par
    #   x      <- 10^tt[(imax.new+1):length(tt)]
    # } else {
    #   loc    <- T.max + loc.par
    #   x      <- tt[(imax.new+1):length(tt)]
    # }
    
    x      <- tt[(imax.new+1):length(tt)]
    
    y.x      <- 1/(1+exp((x - loc)*scale)) + (1 - 1/(1+exp((x[1]-loc)*scale)))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  # five <- wgts.obs[,reconst_fcts[i]]
  # ten <- wgts.obs[,reconst_fcts[i]]
  # fifteen <- wgts.obs[,reconst_fcts[i]]
  # twenty <- wgts.obs[,reconst_fcts[i]]
  # inf <- wgts.obs[,reconst_fcts[i]]
  # wgts0 <- wgts.obs[,reconst_fcts[i]]
  
  # temp <- which(tt %in% t.points)
  # pdf(file = paste0("useful-pics/weight.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(tt, wgts.obs[,reconst_fcts[1]], type='l', lwd=3, col='darkorange',
  #      main="(b) Weight", ylim=c(0,1), xlab="t", ylab="w(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(tt[temp], wgts.obs[temp,reconst_fcts[1]],pch=16, col='darkorange')
  # grid()
  # dev.off()
  
  # save(five, ten, fifteen, twenty, inf, wgts0, file='useful-pics/info_plot_wgts.RData')
  # load('Results/images/info_plot_wgts.RData')
  # load('Results/images/plot-recwgts-info.RData')
  # library(ggplot2)
  # library(scales)
  # gc.ramp <- hue_pal()(8)
  # temp <- which(tt %in% t.points)
  # pdf(file = paste0("Results/images/weight-revision.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(tt, ten, type='l', lwd=3, lty=1, col=gc.ramp[3],
  #      main="(b) Weights", ylim=c(0,1), xlab="t", ylab="w(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(tt[temp], ten[temp],pch=16, col=gc.ramp[3])
  # lines(tt, five, lwd=3, lty=2, col=gc.ramp[2])
  # lines(tt, fifteen, lwd=3, lty=3, col=gc.ramp[4])
  # lines(tt, twenty, lwd=3, lty=4, col=gc.ramp[5])
  # lines(tt, inf, lwd=3, lty=5, col=gc.ramp[6])
  # lines(tt, wgts0, lwd=3, lty=6, col=gc.ramp[8])
  # points(t.points, wgts.obs.KLAL, pch=17, col=gc.ramp[7])
  # legend(0,0.75, legend=c('a=5', 'a=10', 'a=15', 'a=20', 'a=inf', '0-wgts','rec-wgts'),
  #        col=c(gc.ramp[2],gc.ramp[3],gc.ramp[4],gc.ramp[5],gc.ramp[6],gc.ramp[8],gc.ramp[7]),
  #        lty=c(1,2,3,4,5,6,NA),pch=c(NA,NA,NA,NA,NA,NA,17), lwd=c(3,3,3,3,3,3,NA), cex=1.6)
  # grid()
  # dev.off()
  
  # 
  ## Weights as functional data
  ## Smooth the weights on a B-spline basis of order 2
  basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=2)
  wgts.fd  <- smooth.basis(tt, wgts.obs, basis)$fd
  
  ## Output
  out        <- list(wgts.obs, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

#### Definition of zero weights
create_zero_weights <- function(curves.rec, t.points, breaks, reconst_fcts, Thp, log.flag=FALSE)
{
  ## Utilities
  n <- dim(curves.rec)[2]
  
  ## Building the weights
  if(log.flag)
  {
    tt <- t.points
  } else {
    step <- round(min(t.points[2:length(t.points)] - t.points[1:(length(t.points)-1)]), digits=2)
    tt <- seq(range(t.points)[1], range(t.points)[2], by=step)
  }
  
  
  wgts.obs <- matrix(data=1, nrow=length(tt), ncol=n)
  
  for(i in 1:length(reconst_fcts))
  {
    t.max    <- Thp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    
    x        <- tt[(imax.new+1):length(tt)]
    
    y.x      <- rep(1e-6, length(x))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  ## Weights as functional data
  ## Smooth the weights on a B-spline basis of degree 2
  basis   <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=1)
  wgts.fd <- smooth.basis(tt, wgts.obs, basis)$fd
  
  ## Output
  
  out        <- list(wgts.obs, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

## Prediction
predict_fRegress <- function(mod, xlist, t.points)
{
  
  ## Utilities
  y         <- mod$yfdobj
  p         <- length(xlist)
  n         <- dim(xlist[[1]]$coefs)[2]
  rangeval  <- range(t.points)
  ybasisobj <- y$basis
  ynbasis   <- ybasisobj$nbasis
  nfine     <- max(501,10*ynbasis+1)
  tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
  
  ## Evaluate it component per component and function by function
  yhatmat  <- matrix(0, nrow=nfine, ncol=n)
  for (j in 1:p) {
    xfdj       <- eval.fd(tfine, xlist[[j]], 0, FALSE)
    betafdParj <- mod$betaestlist[[j]]
    betafdj    <- betafdParj$fd
    betavecj   <- eval.fd(tfine, betafdj, 0, FALSE)
    yhatmat    <- yhatmat + as.vector(betavecj)*xfdj
  }
  yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
  
  ## Output
  return(yhatfdobj)
  
}

## Error evaluation -----------------------------------------------------------
eval_MSE_functional <- function(curves, curves.hat, t.points)
{
  n <- dim(curves$coefs)[2]
  
  delta.T <- range(t.points)[2] - range(t.points)[1]
  MSE.vec <- numeric(n)
  for(i in 1:n) #i=1
  {
    yi <- curves
    yi$coefs <- as.matrix(curves$coefs[,i], ncol=1)
    
    yi.hat <-  curves.hat
    yi.hat$coefs <- as.matrix(curves.hat$coefs[,i], ncol=1)
    
    err.i <- yi - yi.hat
    MSE.vec[i] <- diag(inprod(err.i, err.i, 0, 0, rng=range(t.points)))
    MSE.vec[i] <- MSE.vec[i]/delta.T
  }
  
  MSE     <- sum(MSE.vec)/n
  
  return(MSE)
}

eval_MSE <- function(curves, curves.hat, t.points)
{
  ## Utilities
  n <- dim(curves)[2]
  
  ## MSE_b
  ## per ogni funzione nel test
  ## valuto la distanza di quella predetta da quella vera, poi integro
  ## ma soltanto nella parte di dominio dove questa funzione è osservata
  MSE <- 0
  MSE_reconstruction <- 0
  
  range_reconstruction <- c(log10(5),log10(10))
  
  n.complete <- 0
  
  for(j in 1:n)
  {
    curve       <- curves[,j]
    M_bool_vec  <- is.na(curve)
    O_bool_vec  <- !M_bool_vec
    
    t.O         <- t.points[O_bool_vec]
    curve.O     <- curve[O_bool_vec]
    
    curve.hat.O <- curves.hat[O_bool_vec,j]
    
    err.O       <- curve.O - curve.hat.O
    
    basis    <- create.bspline.basis(rangeval=range(t.O), breaks=t.O, norder=4)
    fPar     <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    err.O.fd <- smooth.basis(t.O, err.O, fPar)$fd
    
    MSE_j    <- inprod(err.O.fd, err.O.fd, 0, 0, rng=range(t.O))
    
    Tf <- tail(t.O, n=1)
    Ti <- head(t.O, n=1)
    
    delta.T <- Tf - Ti
    MSE     <- MSE + MSE_j/delta.T
    
    if(sum(O_bool_vec)==length(t.points))
    {
      MSE_reconstruction_j <- inprod(err.O.fd, err.O.fd, 0, 0, rng=range_reconstruction)
      MSE_reconstruction <- MSE_reconstruction + MSE_reconstruction_j/(range_reconstruction[2]-range_reconstruction[1])
      n.complete <- n.complete + 1
    }
  }
  
  MSE <- MSE/n
  MSE_reconstruction <- MSE_reconstruction/n.complete
  
  out <- list(MSE                = MSE,
              MSE_reconstruction = MSE_reconstruction)
  
  return(out)
}
