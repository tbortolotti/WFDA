setwd('~/Documents/R/WFDA')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(R.matlab)
library(latex2exp)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(snowfall)
library(psych)
library(progress)
library(beepr)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data ---------------------------------------------------------------
load('output/preprocessed-data/xlist-logg.RData')
load('output/preprocessed-data/data.RData')

## Load functions
source('code/functions/weighted-analysis.R')
source('code/functions/weighted-regression.R')
source('code/functions/weighted-smoothing.R')
source('code/functions/reconstruction.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
q <- length(xlist)
reconst_fcts  <- find_obs_inc(Y = curves)

fix.par <- 100
t.points <- log10(T.period)
t.points[1] <- -2.5
breaks <- t.points

data <- list(dJB  = dJB,
             MAG  = MAG,
             SoF  = SoF,
             VS30 = VS30)

data.f <- data.frame(dJB = dJB, MAG = MAG, SoF = SoF, VS30 = VS30)

(Start.Time <- Sys.time())
## Extrapolation ---------------------------------------------------------------
extrapolate   <- extrapolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.extrap <- extrapolate$curves.rec

## Construction of the weights -------------------------------------------------
wgt       <- create_weights(curves.rec    = curves.extrap,
                            t.points      = t.points,
                            breaks        = breaks,
                            fix.par       = fix.par,
                            reconst_fcts  = reconst_fcts,
                            Thp           = log10(T_hp))

## Smoothing -------------------------------------------------------------------
smth             <- wt_bsplinesmoothing(curves   = curves.extrap,
                                        wgts.obs = wgt$wgts.obs,
                                        t.points = t.points,
                                        breaks   = breaks,
                                        lambda   = 1e-5,
                                        set.cb   = FALSE)
curves.extrap.fd <- smth$curves.fd

## B-list ----------------------------------------------------------------------
load('output/calibration/blist_EAASS_multistart.RData')

## Regression and beta estimation ----------------------------------------------
mod <- weighted_fRegress(y            = curves.extrap.fd,
                         xfdlist      = xlist,
                         betalist     = blist,
                         wgts         = wgt$wgts.fd)

## y_hat   ---------------------------------------------------------------------
curves.extrap.hat   <- predict_fRegress(mod          = mod,
                                        xlist        = mod$xfdlist,
                                        t.points     = t.points)
curves.extrap.hat.v <- eval.fd(t.points, curves.extrap.hat)
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

## MSE evaluation and Sigma estimation -----------------------------------------
#' through an event-wise crossvalidation
pw_MSE_Sigma.val <- pwMSE(curves    = curves,
                          curves.fd = curves.extrap.fd,
                          xlist     = xlist,
                          t.points  = t.points,
                          events    = event.id,
                          blist     = blist,
                          B         = 10,
                          wgts.fd   = wgt$wgts.fd,
                          set.ITA18 = TRUE,
                          data      = data.f,
                          wgts.flag = TRUE)

beep()
save(pw_MSE_Sigma.val, file='output/application/pw_MSE_Sigma.RData')
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
beep()

## Model comparison ------------------------------------------------------------

source('code/functions/application-plots.R')
source('code/functions/scalar-analysis.R')

## Fit scalar ITA18
scalar.ITA18 <- fit_ITA18(data.f, curves)

coefs.ITA18 <- scalar.ITA18$coefs.ITA18
y.hat.ITA18 <- scalar.ITA18$y.hat.ITA18

save(coefs.ITA18, file="output/application/coefs_ITA18.RData")

## Plot MSE and sigma comparison
load('output/application/pw_MSE_Sigma.RData')

name_dir     <- paste0("output/images/model-comparison")

plotMSE_pw(MSE.vec   = pw_MSE_Sigma.val$MSE_t,
           MSE.ita18 = pw_MSE_Sigma.val$MSE_t18,
           t.points  = t.points,
           name_dir  = name_dir)

plot_sigma(sigma.t     = pw_MSE_Sigma.val$Sigma_t,
           sigma.ITA18 = pw_MSE_Sigma.val$Sigma_t18,
           t.points    = t.points,
           name_dir    = name_dir)


## Bootstrap sample of the functional coefficients -----------------------------
## Preprocessing
extrapolate   <- extrapolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.extrap <- extrapolate$curves.rec

## Construction of the weights
wgt       <- create_weights(curves.rec    = curves.extrap,
                            t.points      = t.points,
                            breaks        = breaks,
                            fix.par       = fix.par,
                            reconst_fcts  = reconst_fcts,
                            Thp           = log10(T_hp))

## Smoothing
smth             <- wt_bsplinesmoothing(curves   = curves.extrap,
                                        wgts.obs = wgt$wgts.obs,
                                        t.points = t.points,
                                        breaks   = breaks,
                                        lambda   = 1e-5,
                                        set.cb   = FALSE)
curves.extrap.fd <- smth$curves.fd
L <- curves.extrap.fd$basis$nbasis

## Generate a bootstrap sample of functional coefficients ----------------------

# 1. Fit the regression
mod <- weighted_fRegress(y            = curves.extrap.fd,
                         xfdlist      = xlist,
                         betalist     = blist,
                         wgts         = wgt$wgts.fd)

# 2. evaluate the residuals
curves.extrap.hat <- predict_fRegress(mod          = mod,
                                      xlist        = mod$xfdlist,
                                      t.points     = t.points)
res <- curves.extrap.hat - curves.extrap.fd
wgts.fd <- wgt$wgts.fd

# 3. Repeatedly sample from the empirical distribution of data
set.seed(140996)
B <- 1000   # Time consuming. One may decide to lower B and do some trials

obs <- seq(1,n)
B.list <- array(data=0, dim=c(N,q,B))

B.mat <- matrix(nrow=N, ncol=q)

for(b in 1:B)
{
  print(paste0("Iteration ", b, " of ", B))
  obs.b <- sample(obs, replace=T)
  
  res.b <- res
  res.b$coefs <- res$coefs[,obs.b]
  curves.extrap.fd.b <- curves.extrap.hat + res.b
  
  mod.b <- weighted_fRegress(y            = curves.extrap.fd.b,
                             xfdlist      = xlist,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
  
  for(i in 1:q)
  {
    B.mat[,i] <- eval.fd(t.points, mod.b$betaestlist[[i]]$fd)
  }
  
  B.list[,,b] <- B.mat
  
}

beep()

save(B.list, mod, n, q, L, N, T.period, file='output/application/bootstrap_result_B1000.RData')


## PLOTS ------------------------------------------------------------------------------
load("output/application/coefs_ITA18.RData")
a0      <- coefs.ITA18$a0
b1      <- coefs.ITA18$b1
b2      <- coefs.ITA18$b2
c1      <- coefs.ITA18$c1
c2      <- coefs.ITA18$c2
c3      <- coefs.ITA18$c3
k0      <- coefs.ITA18$k0
f1      <- coefs.ITA18$f1
f2      <- coefs.ITA18$f2
X.ita18 <- cbind(a0,b1,b2,f1,f2,c1,c2,c3,k0)

load('output/application/bootstrap_result_B1000.RData')
library(boot)

overall.conf <- list()
for(j in 1:q)
{
  overall.conf[[j]] <- envelope(mat = t(B.list[,j,]), level=0.95)$overall
}

## Regression coefficients boxplots and ITA18 comparison
library(wesanderson)

pal <- wes_palette('Cavalcanti1')
rect.col1 <- adjustcolor(pal[2], alpha=0.1)
pal.gb <- "grey70"
pal.gb.fade <- adjustcolor(pal.gb, alpha = 0.5)



t.points <- log10(T.period)
t.points[1] <- -2.5
xtick <- seq(-2,1,by=1)

names.coefs <- c('a','b1','b2','f1','f2','c1','c2','c3','k')
coefs.names <- c('Intercept', 'B1', 'B2', 'F1', 'F2','C1','C2','C3','K')

for(j in 1:q) #j=2
{
  {
    beta.j.mod <- mod$betaestlist[[j]]$fd
    beta.j.val <- eval.fd(t.points, beta.j.mod)
    
    y.min <- min(c(overall.conf[[j]][2,], X.ita18[,j]))
    y.max <- max(c(overall.conf[[j]][1,], X.ita18[,j]))
    ylim <- c(y.min, y.max)
    
    pdf(file = paste0("output/images/model-comparison/reg_",coefs.names[j],".pdf"), width = 8, height = 5)
    par(mar=c(4.5, 4, 2.5, 1)+.1)
    fda::fbplot(fit=B.list[,j,], x=t.points, method="MBD", color=pal.gb.fade,barcol="grey60", xlab='Period [s]',
                xlim=range(t.points), ylim=ylim, cex.axis=1.8, cex.lab=1.8, ylab='', xaxt='n')
    title(main=paste0('(a) ',names.coefs[j]), cex.main=1.8)
    points(t.points, X.ita18[,j], pch=16, col="black")
    grid()
    abline(h=0, col="black", lty=4, lwd=1)
    axis(side=1, at=xtick, labels = 10^xtick, cex.axis=2)
    rect(xleft=t.points[21],
         ybottom=ylim[1],
         xright=t.points[37],
         ytop=ylim[2],
         col=rect.col1,
         border=NA)
    dev.off()
  }
  
}

