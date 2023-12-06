library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(psych)
library(progress)
library(fields)
library(ggplot2)
library(beepr)
library(matrixcalc)
library(latex2exp)
library(ggplotify)

rm(list=ls())
graphics.off()
cat("\014")

#' We recommend users to take advantage of the Outline tool available in R, in
#' order to move easily in the script.

## Load Functions
source('code/functions/weighted-analysis.R')
source('code/functions/generate-data.R')
load('output/simulation/reg_info.RData')

## Simulation - SET 2  --------------------------------------------------------------
perc.po <- 4
fix.par <- "unwgt"

## Set the case information
case.info <- list(n.sim       = 101,
                  ext.noise   = 0.5,
                  perc        = 0.1*perc.po,
                  left.bound  = 1.5,
                  right.bound = 3.5)

B <- 100
MSE <- numeric(B)

method <- 'KLAl'

#### Unweighted -------------------------------------------------------------

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = FALSE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'_nowgts.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### a=5 -------------------------------------------------------------

fix.par <- 5

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-par', fix.par, '.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### a=10 -------------------------------------------------------------

fix.par <- 10

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-par', fix.par, '.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### a=15 -------------------------------------------------------------

fix.par <- 15

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-par', fix.par, '.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

#### a=20 -------------------------------------------------------------

fix.par <- 20

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-par', fix.par, '.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### a=inf -------------------------------------------------------------

fix.par <- 100

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-par', fix.par, '.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### Correlation-driven weights -------------------------------------------------------------

fix.par <- 100

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = TRUE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-corwgts.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### 0 weights -------------------------------------------------------------

fix.par <- "0-wgts"

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  curves.true <- simulated_data$curves.true
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b               = 1,
                                                  B               = case.info$n.sim,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

name.file <- paste0('output/simulation/',method,'-0wgts.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

## Plot results -------------------------------------------------------
#' With the following chunks of code, we are replicating Figure 5 of the
#' manuscript (a and b). Specifically, Subsection MSE of this script reproduces
#' Figure 5a, while Subsection Variance reproduces Figure 5b

beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3 
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

box.dir <- 'output/images/simulation'


#### SET 2 - MSE ----------------------------------------------------
levs <- 8
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights","a=5","a=10","a=15", "a=20", "a=inf", "rec-wgts", "0-wgts"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

{
  load('output/simulation/KLAl_nowgts.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par5.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par15.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par20.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par100.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  
  load('output/simulation/KLAl-corwgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  
  load('output/simulation/KLAl-0wgts.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
}


# Defining the dataframe for the plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights","a=5","a=10","a=15" ,"a=20", "a=inf", "rec-wgts", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() + scale_fill_grey()
pgplot <- pgplot +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() + 
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'), fill = "Parameter: ", title="(a) Weights definition: MSE") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
pgplot + 
  theme(legend.position="none")
ggsave(filename = "set2-mse.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)

pg_legend <- cowplot::get_legend(pgplot)
as.ggplot(pg_legend)


ggsave(filename = "set2-legend.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 14,
       height = 0.7,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)


#### SET 2 - Variance ------------------------------------------------
levs <- 8
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=5", "a=10", "a=15", "a=20", "a=inf", "rec-wgts", "0-wgts"),
                   each=B)
beta_var   <- numeric(3*levs*B)

{
  load('output/simulation/KLAl_nowgts.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par5.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par15.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par20.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-par100.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  
  load('output/simulation/KLAl-corwgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  load('output/simulation/KLAl-0wgts.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
}

# Defining the dataframe for the plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("No weights", "a=5", "a=10", "a=15", "a=20", "a=inf", "rec-wgts", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

## Grouped boxplot
pg_plot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "Parameter", title="(b) Weights definition: Variance")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))
pg_plot + theme(legend.position = "none")

ggsave(filename = "set2-var.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## LOO CROSS-VALIDATION ---------------------------------------------------------
# Simulate the data
# Set the case information
case.info <- list(n.sim       = 100,
                  ext.noise   = 0.5,
                  perc        = 0.4,
                  left.bound  = 1.5,
                  right.bound = 3.5)

seed <- 140996
simulated_data <- generate_data(seed      = seed,
                                reg.info  = reg.info,
                                case.info = case.info)

t.points       <- simulated_data$t.points
T_hp           <- simulated_data$T_hp
curves         <- simulated_data$curves
curves.true.fd <- simulated_data$curves.true.fd
xlist          <- simulated_data$xlist

## Utilities
B <- 100
method <- 'KLAl'

#### Unweighted -------------------------------------------------------------
fix.par <- 0
MSE_cv <- matrix(data=0, nrow=B, ncol=1)

(Start.Time <- Sys.time())

print(paste0('Parameter ', fix.par))
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  method_evaluation <- workflow_weighted_analysis(b               = b,
                                                  B               = B,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = FALSE,
                                                  wgts.recon.flag = FALSE)
  
  MSE_cv[b,1] <- method_evaluation$MSE
  pb$tick()
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/simulation/LOO_',method,'_nowgts.RData')
save(MSE_cv, file=name.file)


#### Logistic weights ------------------------------------------------------- 

vec.par <- list(5,10,15,20,100,'0-wgts')
MSE_cv <- matrix(data=0, nrow=B, ncol=length(vec.par))

(Start.Time <- Sys.time())
for(a in 1:length(vec.par)){ # a <- 6
  fix.par <- vec.par[[a]]
  print(paste0('Parameter ', fix.par))
  pb <- progress_bar$new(total=B)
  for(b in 1:B) # b=1
  {
    method_evaluation <- workflow_weighted_analysis(b               = b,
                                                    B               = B,
                                                    t.points        = t.points,
                                                    breaks          = t.points,
                                                    T_hp            = T_hp,
                                                    curves          = curves,
                                                    curves.true.fd  = curves.true.fd,
                                                    xlist           = xlist,
                                                    method          = method,
                                                    fix.par         = fix.par,
                                                    wgts.flag       = TRUE,
                                                    wgts.recon.flag = FALSE)
    
    MSE_cv[b,a] <- method_evaluation$MSE
    pb$tick()
  }
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/simulation/LOO_',method,'_wgts.RData')
save(MSE_cv, file=name.file)

#### Reconstruction-driven weights ----------------------------------------------
fix.par <- 0
MSE_cv <- matrix(data=0, nrow=B, ncol=1)

(Start.Time <- Sys.time())

pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  method_evaluation <- workflow_weighted_analysis(b               = b,
                                                  B               = B,
                                                  t.points        = t.points,
                                                  breaks          = t.points,
                                                  T_hp            = T_hp,
                                                  curves          = curves,
                                                  curves.true.fd  = curves.true.fd,
                                                  xlist           = xlist,
                                                  method          = method,
                                                  fix.par         = fix.par,
                                                  wgts.flag       = TRUE,
                                                  wgts.recon.flag = TRUE)
  
  MSE_cv[b,1] <- method_evaluation$MSE
  pb$tick()
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('output/simulation/LOO_',method,'_corwgts.RData')
save(MSE_cv, file=name.file)

# Plot results-------------------------------------------------------------
#' With the following chunks of code, we are replicating Figure 6 of the
#' manuscript.

box.dir <- 'output/simulation/images'

MSE_cv.mat <- matrix(data=0, nrow=B, ncol=1)
load('output/simulation/LOO_KLAl_nowgts.RData')
MSE_cv.mat[,1] <- MSE_cv

load('output/simulation/LOO_KLAl_wgts.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,1:5])

load('output/simulation/LOO_KLAl_corwgts.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,1])

load('output/simulation/LOO_KLAl_wgts.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,6])

MSE <- vec(MSE_cv.mat)

class <- c(rep("unwgt", B),
           rep('a=5', B),
           rep('a=10', B),
           rep('a=15', B),
           rep('a=20', B),
           rep('a=inf', B),
           rep("rec-wgts",B),
           rep('0-wgts', B))

data.box <- data.frame(MSE, as.factor(class))
names(data.box) <- c('MSE', 'class')
class_order<- c("unwgt", "a=5","a=10", "a=15", "a=20", "a=inf", "rec-wgts", "0-wgts")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

box.dir <- 'output/images/simulation'

pg_plot <- ggplot(data.box, aes(x = class, y = MSE, fill=class))+
  geom_boxplot() + scale_fill_grey()
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.04)) +
  theme(legend.position="None") + 
  theme_bw() +
  labs(x="", y=TeX(r'($ || \hat{y}_{i} - y_i ||_2^2$)'), title="Weights definition: LOO MSE") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

ggsave(filename = "set2-LOOmse.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 10,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)

dev.off()
