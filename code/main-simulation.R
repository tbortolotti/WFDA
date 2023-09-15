setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

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

rm(list=ls())
graphics.off()
cat("\014")

## Load Functions --------------------------------------------------------------
source('Simulation/methods/workflow_weighted_analysis_rev.R')
source('Simulation/methods/generate_data.R')
load('Simulation/DATA/reg_info.RData')

## Simulation ------------------------------------------------------------------
fix.par <- 10
perc.po <- 4
smooth.noise <- 10

## Set the case information
case.info <- list(n.sim       = 101,
                  noise       = 0.01*smooth.noise, #all analysis with noise=0.1 (i.e. smooth.noise=10)
                  ext.noise   = 0.5,
                  perc        = 0.1*perc.po,
                  left.bound  = 1.5,
                  right.bound = 3.5)

B <- 100
MSE<- numeric(B)

method <- 'Kraus'
#method <- 'KLNoAl'
#method <- 'KLAl'

#b=1
(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  case      = "CASE-1",
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

name.file <- paste0('Simulation/Results/repeated-simulations/',method,'-corwgts-new.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

