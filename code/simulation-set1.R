
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

## Simulation - SET 1 --------------------------------------------------------
fix.par <- 10
perc.po <- 4

## Set the case information
case.info <- list(n.sim       = 101,
                  ext.noise   = 0.5,
                  perc        = 0.1*perc.po,
                  left.bound  = 1.5,
                  right.bound = 3.5)

B <- 100
MSE<- numeric(B)

#### Kraus, weighted analysis ------------------------------------------

method <- 'Kraus'

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

name.file <- paste0('output/simulation/',method,'.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

#### KL-PC, weighted analysis -------------------------------------------------
method <- 'KLNoAl'

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

name.file <- paste0('output/simulation/',method,'.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

#### KL-AL, weighted analysis -------------------------------------------------
method <- 'KLAl'

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

name.file <- paste0('output/simulation/',method,'-par10.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)


#### Kraus, unweighted analysis ------------------------------------------

method <- 'Kraus'

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

#### KL-PC, unweighted analysis -------------------------------------------------
method <- 'KLNoAl'

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

#### KL-AL, unweighted analysis -------------------------------------------------
method <- 'KLAl'

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


## Plot results -------------------------------------------------------
#' With the followig chunks of code, we are replicating Figure 4 of the
#' manuscript (a and b), and the values of bias^2 entering Table 1 of the
#' manuscript. Specifically, Subsection MSE of this script reproduces
#' Figure 4a, Subsection Variance reproduces Figure 4b, and Subsection bias^2
#' reproduces the values in Table 1

beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3 
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

box.dir <- 'output/images/simulation'


#### MSE ---------------------------------
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

{
  load('output/simulation/Kraus.RData')
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
  
  load('output/simulation/KLNoAl.RData')
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
  
  load('output/simulation/Kraus_nowgts.RData')
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
  
  load('output/simulation/KLNoAl_nowgts.RData')
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
  
  load('output/simulation/KLAl_nowgts.RData')
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
  
}

# Defining the dataframe for the plot
data <- data.frame(coefficient, method, beta_mse)
method_order  <- c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'),
       fill = "Method: ", title="(a) Reconstruction methods: MSE")+
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
ggsave(filename = "set1-MSE.pdf",
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
as.ggplot(pg_legend)

ggsave(filename = "set1-legend.pdf",
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
dev.off()

#### Variance -----------------------
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("Kraus: wgt", "KL-PC: wgt", "KL-AL: wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=B)
beta_var   <- numeric(3*levs*B)
bias0.vec   <- numeric(levs)
bias1.vec   <- numeric(levs)
bias2.vec   <- numeric(levs)

{
  load('output/simulation/Kraus.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLNoAl.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
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
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/Kraus_nowgts.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLNoAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLAl_nowgts.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
}

# Defining the dataframe for the plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("Kraus: wgt", "KL-PC: wgt", "KL-AL: wgt", "Kraus", "KL-PC", "KL-AL")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() + scale_fill_grey() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "Method", title="(b) Reconstruction methods: Variance")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))
pgplot + theme(legend.position = "none")

ggsave(filename = "set1-var.pdf",
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

#### bias^2 ----------------------------------------------------------------------------
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs)
method      <- rep(c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=1)
bias2   <- numeric(3*levs)

{
  load('output/simulation/Kraus.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLNoAl.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/Kraus_nowgts.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLNoAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('output/simulation/KLAl_nowgts.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  bias2
  
}
