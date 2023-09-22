setwd('~/Documents/R/WFDA')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(latex2exp)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(snowfall)
library(psych)
library(progress)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data
load('output/preprocessed-data/xlist-logg.RData')
load('output/preprocessed-data/data.RData')

## Load functions
source('code/functions/weighted-analysis.R')
source('code/functions/weighted-regression.R')
source('code/functions/reconstruction.R')

## Utilities
n            <- dim(curves)[2]
q            <- length(xlist)
reconst_fcts <- find_obs_inc(Y = curves)
t.points     <- log10(T.period)
t.points[1]  <- -2.5
breaks       <- t.points

## Selection of the penalization parameters ------------------------------------
## Maintain only the fully observed curves
curves.full <- curves[,-reconst_fcts]
event.id    <- event.id[-reconst_fcts]
xlist.full  <- list()

for(i in 1:length(xlist))
{
  if (inherits(xlist[[i]], "fd"))
  {
    xlist.full[[i]]  <- xlist[[i]]
    xlist.full[[i]]$coefs  <- xlist[[i]]$coefs[,-reconst_fcts]
    
  } else if (inherits(xlist[[i]], "numeric")) {
    
    xlist.full[[i]]  <- xlist[[i]][-reconst_fcts]
    
  } else if (inherits(xlist[[i]], "matrix" )) {
    
    xlist.full[[i]]  <- xlist[[i]][-reconst_fcts,1]
    
  }
}

## Unweighted Smoothing
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-7,1, by=1)
lambda.vec <- sort(10^esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, curves.full, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
curves.fd <- smooth.basis(t.points, curves.full, fPar)$fd

#### Multistart EAASS algorithm -------------------------------------------


calibrate <- FALSE

if(calibrate){
  # Algorithm parameters
  n.iter <- 15
  B <- 5
  P.size <- 20
  n.par <- length(xlist.full)
  r=0.5
  blist.default <- list(fPar, fPar, fPar, fPar, fPar, fPar, fPar, fPar, fPar)
  perturbation.vec <- c(0.5,2)
  
  n.multistart <- 5
  seed.in <- 22
  
  for(b in 1:n.multistart){ #b <- 1
    seed <- seed.in + b
    print(paste0("Iteration ",b, " of ", n.multistart))
    
    set.seed(seed)
    # Initial population of combinations of tuning parameters
    exp.range <- seq(-5,2,by=1)
    lambda.range <- 10^exp.range
    P.initial <- matrix(data=NA, nrow=P.size, ncol=n.par)
    
    for(i in 1:P.size){
      P.initial[i,] <- sample(lambda.range, size=n.par, replace=TRUE)
    }
    
    # Compute the errors for the initial P via a B-fold CV
    # Errors corresponding to each row of P.initial are saved in vector V.initial
    (Start.Time <- Sys.time())
    V.initial <- c()
    pb <- progress_bar$new(total=(P.size), format = "  computing [:bar] :percent eta: :eta")
    for(i in 1:P.size){
      # define blist
      blist <- blist.default
      for(reg in 1:n.par){
        blist[[reg]]$lambda <- P.initial[i,reg] 
      }
      
      # evaluate the MSE via B-fold CV
      MSE_pw <- pwMSE(curves    = curves.full,
                      curves.fd = curves.fd,
                      xlist     = xlist.full,
                      t.points  = t.points,
                      events    = event.id,
                      blist     = blist,
                      B         = 10,
                      wgts.fd   = NULL,
                      wgts.flag = FALSE,
                      verbose   = FALSE)$MSE_t
      V.initial[i] <- sum(MSE_pw^2)/length(MSE_pw)
      
      pb$tick()
    }
    
    End.Time <- Sys.time()
    round(End.Time - Start.Time, 2)
    #beep()
    
    save(P.initial, V.initial, file=paste0("data/calibration/Initial_blist_select_seed",seed,".RData"))
    #load(paste0("data/calibration/Initial_blist_select_seed",seed,".RData"))
    
    V <- V.initial
    P <- P.initial
    
    #set.seed(b)
    (Start.Time <- Sys.time())
    for(k in 1:n.iter){
      print(paste0("Iteration: ",k))
      
      V.sort <- sort(V, decreasing=TRUE)
      idxs.in <- which(V %in% V.sort[1:(r*P.size)])
      if(length(idxs.in)>(r*P.size)){ #solve ties by randomly selecting only r*P.size indexes
        idxs.in <- sample(idxs.in, size=(r*P.size), replace=FALSE)
      }
      Q <- P[idxs.in,]
      Z <- V[idxs.in]
      
      perturbation <- sample(perturbation.vec, size=(r*P.size*n.par), replace=TRUE)
      pert_mat <- matrix(perturbation, nrow=(r*P.size), ncol=n.par)
      Q_prime <- Q * pert_mat
    
      # ---
      temp <- Q_prime
      idxs.ext <- which(Q_prime < 1e-5 | Q_prime > 1e2)
      temp[idxs.ext] <- sample(lambda.range, size=length(idxs.ext), replace=TRUE)
      Q_prime <- temp
      # ---
      
      pb <- progress_bar$new(total=(r*P.size), format = "  computing [:bar] :percent eta: :eta")
      Z_prime <- c()
      for(i in 1:(r*P.size)){
        # define blist
        blist <- blist.default
        for(reg in 1:n.par){
          blist[[reg]]$lambda <- Q_prime[i,reg] 
        }
        
        # evaluate the MSE via B-fold CV
        MSE_pw <- pwMSE(curves    = curves.full,
                        curves.fd = curves.fd,
                        xlist     = xlist.full,
                        t.points  = t.points,
                        events    = event.id,
                        blist     = blist,
                        B         = 10,
                        wgts.fd   = NULL,
                        wgts.flag = FALSE,
                        verbose   = FALSE)$MSE_t
        Z_prime[i] <- sum(MSE_pw^2)/length(MSE_pw)
        pb$tick()
      }
      
      P.in <- P[idxs.in,]
      P <- rbind(P.in, Q_prime)
      V.in <- V[idxs.in]
      V <- c(V.in, Z_prime)
    }
    End.Time <- Sys.time()
    round(End.Time - Start.Time, 2)
    #beep()
    
    save(P, V, file=paste0("data/calibration/P_V_blist_select_seed",seed,".RData"))

    best <- which(V==min(V))
    lambda.opt <- P[best,]
    blist <- blist.default
    for(reg in 1:n.par){
      blist[[reg]]$lambda <- lambda.opt[reg]
    }
    
    save(blist, file=paste0('data/calibration/blist_EAASS_seed',seed,'.RData'))
  }
} else {
  ## Load results and define blist
  n.multistart <- 5
  n.par <- length(xlist.full)
  P.full <- matrix(data=NA, nrow=0, ncol=n.par)
  V.full <- c()
  seed.in <- 22
  for(b in 1:n.multistart){ # b <- 1
    seed <- seed.in + b
    load(paste0("output/calibration/P_V_blist_select_seed",seed, ".RData"))
    
    P.full <- rbind(P.full, P)
    V.full <- c(V.full, V)
  }
  
  plot(1:length(V.full), V.full, type='l')
  
  best <- which(V.full==min(V.full))
  lambda.opt <- P.full[best,]
  lambda.opt
  blist <- blist.default
  for(reg in 1:n.par){
    blist[[reg]]$lambda <- lambda.opt[reg]
  }
  
  save(blist, file=paste0('output/calibration/blist_EAASS_multistart.RData'))
  
}


## Function for the evaluation of the MSE for each proposed method
source('methods/Reconstruction/methods_workflow.R')

# Logarithm of the period
t.points    <- log10(T.period)
t.points[1] <- -2.5
N <- length(t.points)
log.Thp     <- log10(T_hp)

## Selection of the weights ----------------------------------------------------
#' Via event-wise cross-validation
B     <- 10
method <- 'KLAl'

#### Unweighted analysis -------------------------------------------------------
vec.par <- c(0)
A <- length(vec.par)

MSE_cv <- array(data=0, dim=c(N, B, A))
MSE_glob.list <- list(list())

(Start.Time <- Sys.time())
for(a in 1:A) # a=1
{
  fix.par <- vec.par[a]
  for(b in 1:B) # b=1
  {
    print(paste0("Par:", a, " b:",b))
    method_evaluation <- weighted_analysis_event(b          = b,
                                                 T.period   = T.period,
                                                 t.points   = t.points,
                                                 breaks     = t.points,
                                                 T_hp       = T_hp,
                                                 log.Thp    = log.Thp,
                                                 curves     = curves,
                                                 events     = event.id,
                                                 B          = B,
                                                 xlist      = xlist,
                                                 blist      = blist,
                                                 method     = method,
                                                 fix.par    = fix.par,
                                                 wgts.flag  = TRUE,
                                                 seed       = 14996)
    
    MSE_cv[,b,a] <- method_evaluation$MSE_pw
    MSE_glob.list[[a]][[b]] <- method_evaluation$MSE_glob
    
  }

}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/calibration/MSE_nowgts.RData')
save(vec.par, MSE_cv, MSE_glob.list, file=name.file)

#### Logistic weights ----------------------------------------------------------
vec.par <- list(5,10,15,20,100,"0-weights")
A <- length(vec.par)

MSE_cv <- array(data=0, dim=c(N, B, A))
MSE_glob.list <- list(list(),list(),list(),list(),list(), list())

(Start.Time <- Sys.time())
for(a in 1:A) # a=1
{
  fix.par <- vec.par[a]
  for(b in 1:B) # b=1
  {
    print(paste0("Par:", a, " b:",b))
    method_evaluation <- weighted_analysis_event(b          = b,
                                                 T.period   = T.period,
                                                 t.points   = t.points,
                                                 breaks     = t.points,
                                                 T_hp       = T_hp,
                                                 log.Thp    = log.Thp,
                                                 curves     = curves,
                                                 events     = event.id,
                                                 B          = B,
                                                 xlist      = xlist,
                                                 blist      = blist,
                                                 method     = method,
                                                 fix.par    = fix.par,
                                                 wgts.flag  = TRUE,
                                                 seed       = 14996)
    
    MSE_cv[,b,a] <- method_evaluation$MSE_pw
    MSE_glob.list[[a]][[b]] <- method_evaluation$MSE_glob
    
  }
  
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/calibration/MSE_wgts.RData')
save(vec.par, MSE_cv, MSE_glob.list, file=name.file)

#### Show results --------------------------------------------------------------
load('output/calibration/MSE_nowgts.RData')
m_0 <- round(mean(unlist(lapply(MSE_glob.list[[1]], mean))), digits=4)
s_0 <- round(sd(unlist(lapply(MSE_glob.list[[1]], mean))), digits=4)
print(paste0("unweighted", " MSE = ", m_0, " sd = ", s_0))

load('output/calibration/MSE_wgts.RData')
for(i in 1:(length(vec.par)-1))
{
  m_i <- round(mean(unlist(lapply(MSE_glob.list[[i]], mean))), digits=4)
  s_i <- round(sd(unlist(lapply(MSE_glob.list[[i]], mean))), digits=4)
  print(paste0("a = ", vec.par[i], " MSE = ", m_i, " sd = ", s_i))
}

m_0 <- round(mean(unlist(lapply(MSE_glob.list[[6]], mean))), digits=4)
s_0 <- round(sd(unlist(lapply(MSE_glob.list[[6]], mean))), digits=4)
print(paste0("0-weights", " MSE = ", m_0, " sd = ", s_0))

## Selection of the reconstruction method --------------------------------------
#' Via event-wise cross-validation
B       <- 10
fix.par <- 100

method <- 'extrapolation'

MSE_cv <- array(data=0, dim=c(N, B, 1))
MSE_glob.list <- list(list())

(Start.Time <- Sys.time())
for(b in 1:B) # b=1
{
  print(paste0("b: ",b))
  method_evaluation <- weighted_analysis_event(b          = b,
                                               T.period   = T.period,
                                               t.points   = t.points,
                                               breaks     = t.points,
                                               T_hp       = T_hp,
                                               log.Thp    = log.Thp,
                                               curves     = curves,
                                               events     = event.id,
                                               B          = B,
                                               xlist      = xlist,
                                               blist      = blist,
                                               method     = method,
                                               fix.par    = fix.par,
                                               wgts.flag  = TRUE,
                                               seed       = 14996)
  
  MSE_cv[,b,1] <- method_evaluation$MSE_pw
  MSE_glob.list[[1]][[b]] <- method_evaluation$MSE_glob
  
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/calibration/MSE_extrap.RData')

save(MSE_cv, MSE_glob.list, file=name.file)

#### Show results --------------------------------------------------------------
# Extrapolation
load('output/calibration/MSE_extrap.RData')
m_0 <- round(mean(unlist(lapply(MSE_glob.list[[1]], mean))), digits=5)
s_0 <- round(sd(unlist(lapply(MSE_glob.list[[1]], mean))), digits=5)
print(paste0("extrapolation, MSE ", m_0, " sd ", s_0))

# KLAl, a=100
load('output/calibration/MSE_wgts.RData')
m_i <- round(mean(unlist(lapply(MSE_glob.list[[5]], mean))), digits=5)
s_i <- round(sd(unlist(lapply(MSE_glob.list[[5]], mean))), digits=5)
print(paste0("a = ", vec.par[5], " MSE = ", m_i, " sd = ", s_i))
