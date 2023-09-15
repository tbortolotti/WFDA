setwd('~/Documents/R/WFDA')

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
source('functions/weighted-analysis.R')
source('functions/generate_data.R')
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

name.file <- paste0('output/simulation/',method,'-corwgts-new.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

## Plot simulation results -------------------------------------------------------

beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3 
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

box.dir <- 'output/simulation/images'


#### SET 2 - MSE ----------------------------------------------------
levs <- 8
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights","a=5","a=10","a=15", "a=20", "a=inf", "rec-wgts", "0-wgts"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

{
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par5.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par15.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par20.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par100.RData')
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
  
  
  load('Simulation/Results/repeated-simulations/KLAl-corwgts.RData')
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
  
  
  load('Simulation/Results/repeated-simulations/KLAl-0wgts.RData')
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
as_ggplot(pg_legend)


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
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par5.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par15.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par20.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par100.RData')
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
  
  
  load('Simulation/Results/repeated-simulations/KLAl-corwgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-0wgts.RData')
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

#### SET 3 - MSE --------------------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
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
# Building the dataframe for the plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'),
       fill = "", title="(a) Varying fraction of PO data: MSE") +
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
ggsave(filename = "PO-percentage.pdf",
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

pg_legend <- cowplot::get_legend(pgplot)
as_ggplot(pg_legend)


ggsave(filename = "set3-legend.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 12,
       height = 0.7,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)



#### SET 3 - Variance ----------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
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

# Creating the dataframe for the plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(b) Varying fraction of PO data: Variance") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))

pgplot + theme(legend.position="none")

ggsave(filename = "set3-var.pdf",
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

#### SET 3 - [0,1.75] ----------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
}
# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(a) Variation of the estimator\n in the first half of domain") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))

pgplot + theme(legend.position="none")

ggsave(filename = "0175.pdf",
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


#### SET 3 - [1.75,3] --------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-par10.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
}

# Building 
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(b) Variation of the estimator \n in the second half of domain") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))
pgplot + theme(legend.position="none")

ggsave(filename = "17535.pdf",
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

