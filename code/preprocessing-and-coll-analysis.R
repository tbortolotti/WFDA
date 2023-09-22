setwd('~/Documents/R/WFDA')

rm(list=ls())
graphics.off()
cat("\014")

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(latex2exp)
library(R.matlab)
library(calculus)
library(tidyverse)
library("xtable")
library(psych)
library(progress)
library(wesanderson)
library(ggplot2)
library(ggcorrplot)
library(reshape2)


## PREPROCESSING ---------------------------------------------------------------

#### Load Data -----------------------------------------------------------------
data <- read.csv("data/ITA18_SA_flatfile.csv", header=TRUE, sep=";")

names(data)

#' Joyner-Boore distance
#' NT: Where the information on the Joyner-Boore distance [km] is missing, the value
#'     is set equal to the epicentral distance [km]
dJB <- data$JB_dist
dJB[is.na(dJB)] <- data$epi_dist[is.na(dJB)]
sum(is.na(dJB))

#' Moment Magnitude
MAG <- data$Mw
sum(is.na(MAG))

#' Style of Faulting
#' NT: Data that have undefined style of faulting (fm_type_code=="U") are discarded
SoF <- as.character(data$fm_type_code)
SoF[data$fm_type_code=="U"] <- NA
sum(is.na(SoF))

#' Share-wave velocity V_S30
#' NT: vs30_m_sec contains the values of vs30 directly measured, while vs30_m_sec_WA
#'     contains the values inferred from topography according to Wald and Allen (2007)
#'     correlation. We define VS30 as the vector containing the measured values of vs30
#'     where they are available, and containing the inferred values where direct measurements
#'     are missing.
VS30 <- data$vs30_m_sec
VS30[is.na(data$vs30_m_sec)] <- data$vs30_m_sec_WA[is.na(data$vs30_m_sec)]
sum(is.na(VS30))


#' Profiles of Spectral Acceleration
#' 
#' 


T.period <- c(0,0.010,0.025,0.040,0.050,0.070,0.100,
              0.150,0.200,0.250,0.300,0.350,0.400,0.450,
              0.500,0.600,0.700,0.750,0.800,0.900,1.000,
              1.200,1.400,1.600,1.800,2.000,2.500,3.000,
              3.500,4.000,4.500,5.000,6.000,7.000,8.000,
              9.000,10.000)

# Full matrix of SA observations
rotD50.m <- matrix(NA, nrow=dim(data)[1], ncol=length(T.period))
rotD50.m[,1] <- data$rotD50_pga
for(j in 221:256){
  rotD50.m[,(j-220+1)] <- data[,j]
}


#' High pass filter frequency of the U and the V component [Hz] 
U_hp <- data$U_hp
V_hp <- data$V_hp
sum(is.na(U_hp))
sum(is.na(V_hp))

#' Matrix of the filtered SA observations according to the high-pass corner frequencies
rotD50.f <- rotD50.m
# filtering
for(i in 1:dim(rotD50.m)[1])
{
  f_max <- max(U_hp[i],V_hp[i])
  T_max <- 1/f_max
  for(j in 1:length(T.period))
  {
    if(T.period[j]>T_max)
    {
      rotD50.f[i,j] <- NA
    }
  }
}

#' Event
# load event characteristics
event.id <- data$event_id
sum(is.na(event.id))

uni <- unique(event.id)
aux <- numeric(length(event.id))
for(i in 1:length(uni))
{
  aux[which(event.id==uni[[i]])] <- i
}
event.id <- aux

event.lat <- data$ev_latitude
event.long <- data$ev_longitude
sum(is.na(event.lat))
sum(is.na(event.long))

#' Station
station.id <- data$station_code
uni.stat <- unique(station.id)
length(uni.stat)

## Remove missing values
#' Indexes of completely missing profiles
na.m <- c()
for(i in 1:dim(rotD50.f)[1])
{
  na.m[i] <- sum(is.na(rotD50.f[i,]))
}

idxs.SA.na <- which(na.m==37)
idxs.covs.na <- which(is.na(dJB) | is.na(MAG) | is.na(SoF) | is.na(VS30) | is.na(U_hp) |
                       is.na(V_hp) | is.na(event.id))
idxs.na <- union(idxs.SA.na, idxs.covs.na)

dJB <- dJB[-idxs.na]
MAG <- MAG[-idxs.na]
SoF <- SoF[-idxs.na]
VS30 <- VS30[-idxs.na]
U_hp <- U_hp[-idxs.na]
V_hp <- V_hp[-idxs.na]
event.id <- event.id[-idxs.na]
length(unique(event.id))
event.lat <- event.lat[-idxs.na]
event.long <- event.long[-idxs.na]
rotD50.f <- rotD50.f[-idxs.na,]
rotD50.m <- rotD50.m[-idxs.na,]

#' Incomplete records analysis
F_hp <- c()
T_hp <- c()
for(i in 1:length(U_hp))
{
  F_hp[i] <- max(V_hp[i],U_hp[i])
  T_hp[i] <- min(1/F_hp[i],10)
}

# Percentage of observations that have T_hp<T.period, as function of T.period
prop <- c()
for(t in 1:length(T.period))
{
  prop[t] <- length(which(T_hp>=T.period[t]))/length(T_hp) 
}

# which observations have incomplete domain
na.f <- c()
for(i in 1:dim(rotD50.f)[1])
{
  na.f[i] <- sum(is.na(rotD50.f[i,]))
}
obs.inc <- which(na.f>0)

curves <- log10(t(rotD50.f))

n <- dim(curves)[2]
obs.comp <- (1:n)[-obs.inc]

save(curves, T.period, rotD50.f, dJB, MAG, VS30, SoF,
     event.id, event.lat, event.long,
     U_hp, V_hp, T_hp, obs.comp, obs.inc, file='output/preprocessed-data/data.RData')

#### Regressors construction ---------------------------------------------------
load('data/ITA18_parameters.RData')
load('output/preprocessed-data/data.RData')

library(wesanderson)
pal <- wes_palette('Cavalcanti1')

na.f <- c()
for(i in 1:dim(curves)[2])
{
  na.f[i] <- sum(is.na(curves[,i]))
}
obs.inc <- which(na.f>0)

# col1 <- rep(pal[2],length(dJB))
# col1[obs.inc] <- "darkorange"
dticks <- seq(-1,2, by=1)

auxi <- dJB
idxs <- which(dJB<0.1)
auxi[idxs] <- 0.1

## Plot of magnitude versus distance -- black and white
cols <- rep(pal[1], length(dJB))
cols[which(SoF=="NF")] <- "gray80"
cols[which(SoF=="SS")] <- "gray60"
cols[which(SoF=="TF")] <- "gray40"

pdf(file = paste0("output/images/MAGvsDIST.pdf"), width = 5, height = 5)
par(mar=c(4.5,4.5,2.5,1)+.1)
plot(log10(auxi), MAG, col=cols, xaxt="n", xlab=TeX("$log_{10} (d_{JB})$"), ylab=TeX("$M_w$"),
     main="(b) Magnitude vs JB distance", cex.lab=1.65, cex.axis=1.65, cex.main=1.65, pch=16)
axis(side=1, at=dticks, labels=10^dticks, cex.axis=1.65)
grid()
legend(-1, 5, legend=c("NF", "SS", "TF"), col=c("gray80", "gray60","gray40"),
       pch=16, cex=1.2)
dev.off()

## 1. Build Mh, Mref and h functional objects
names(ITA18.parameters)

par(mfrow=c(3,1))
plot(T.period, ITA18.parameters$Mh.vec, type='l', xlab='Period', ylab='Mh')
plot(T.period, ITA18.parameters$Mref.vec, type='l', xlab='Period', ylab='Mref')
plot(T.period, ITA18.parameters$h.vec, type='l', xlab='Period', ylab='h')

t.points <- log10(T.period)
t.points[1] <- -2.5
xticks <- c(-2,-1,0,1)

par(mfrow=c(3,1))
plot(t.points, ITA18.parameters$Mh.vec, type='l', xlab='log10(T)', ylab='Mh')
plot(t.points, ITA18.parameters$Mref.vec, type='l', xlab='log10(T)', ylab='Mref')
plot(t.points, ITA18.parameters$h.vec, type='l', xlab='log10(T)', ylab='h')
dev.off()

prop <- c()
for(t in 1:length(T.period))
{
  prop[t] <- length(which(T_hp>=T.period[t]))/length(T_hp) 
}

vert.idxs <- c(21, 33)

pdf(file = paste0("output/images/perc-log-period.pdf"), width = 8, height = 5)
par(mar=c(4.5, 4.5, 2.5, 1)+.1)
plot(t.points, prop, ylim=c(0,1), xlab="Period [s]", ylab="Fraction", type='l', lwd=3,
     col='black', xaxt='n', cex.lab=1.8, cex.main=1.8, cex.axis=1.8,
     main = "(b) Period-wise fraction of observed data")
axis(side=1, at=xticks, labels = 10^xticks, cex.axis=1.8)
points(t.points, prop, pch=16, col='black')
grid()
abline(h=0.75, col='black', lty=5, lwd=1)
abline(v=t.points[21], col='black', lty=4, lwd=1)
abline(v=t.points[32], col='black', lty=4, lwd=1)
dev.off()

yticks <- c(-2,0,2)
pdf(file = paste0("output/images/curves.pdf"), width = 8, height = 5)
par(mar=c(4.5, 4.5, 2.5, 1)+.1)
matplot(t.points, curves, type='l', col='grey80', ylab="IM",
        xlab="Period [s]", xaxt='n', yaxt='n',
        cex.lab=1.8, cex.main=1.8, cex.axis=1.8,
        main = "(a) Curves of IM")
axis(side=1, at=xticks, labels = 10^xticks, cex.axis=1.8)
axis(side=2, at=yticks, labels=10^yticks, cex.axis=1.8)
lines(t.points, curves[,194], type='l', lwd=3, col='black')
points(t.points, curves[,194], pch=16, col='black')
grid()
dev.off()

#### Functional covariates on log(T) -------------------------------------------
t.points <- log10(T.period)
t.points[1] <- -2.5
breaks <- t.points
# Mh
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=1)
Mh.fd <- smooth.basis(t.points, ITA18.parameters$Mh.vec, basis)$fd
plot(Mh.fd)
lines(t.points, ITA18.parameters$Mh.vec, type='l', xlab='Period', ylab='Mh', col='red')
dev.off()
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mh.fd <- smooth.basis(t.points, ITA18.parameters$Mh.vec, fPar)$fd


xticks <- c(-2,-1,0,1)
xx <- seq(-2.5, 1, length.out=1000)
Mh.fd.vals <- eval.fd(xx, Mh.fd)
pdf(file = paste0("output/images/Mh.pdf"), width = 8, height = 5)
par(mar=c(4.5, 4.5, 2.5, 1)+.1)
plot(xx, Mh.fd.vals, type='l', xlab="Period [s]", ylab=TeX(r'($M_h$)'),
     lwd=3, col='black', xaxt='n', main='(a) Smoothing of the hinge magnitude',
     cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
axis(side=1, at=xticks, labels = 10^xticks, cex.axis=1.8)
lines(t.points, ITA18.parameters$Mh.vec, type='l', col='grey40', lwd=2, lty=5)
legend(x = -2.5, y=6.2,   # Coordinates (x also accepts keywords)
       legend = c('Lanzano et al., 2019','Quadratic function'),
       col = c("grey40", "black"),
       lwd = 2,
       lty = c(5,1),
       cex = 1.8)
grid()
dev.off()

# Mref
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,10, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.parameters$Mref.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

plot(log10(lambda.vec), gcv.vec, type='l')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
Mref.fd <- smooth.basis(t.points, ITA18.parameters$Mref.vec, fPar)
plot(Mref.fd)
lines(t.points, ITA18.parameters$Mref.vec, type='l', xlab='Period', ylab='Mref', col='red')

# h
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.parameters$h.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

plot(log10(lambda.vec)[2:5], gcv.vec[2:5], type='l')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
h.fd <- smooth.basis(t.points, ITA18.parameters$h.vec, fPar)
plot(h.fd)
lines(t.points, ITA18.parameters$h.vec, type='l', xlab='Period', ylab='h', col='red')
dev.off()

## 2. Create the list of functional regressors
load('output/preprocessed-data/data.RData')

n <- dim(rotD50.f)[1]
N <- length(t.points)

# Create the matrix of values assumed by the regressors at the sampling instants
reg.Mlow <- matrix(data=0, nrow=N, ncol=n)
reg.Mhigh <- matrix(data=0, nrow=N, ncol=n)
reg.D1 <- matrix(data=0, nrow=N, ncol=n)
reg.D2 <- matrix(data=0, nrow=N, ncol=n)
reg.D3 <- matrix(data=0, nrow=N, ncol=n)

R.mat <- matrix(data=0, nrow=N, ncol=n)

for(t in 1:N)
{
  Mh   <- ITA18.parameters$Mh.vec[t]
  Mref <- ITA18.parameters$Mref.vec[t]
  h    <- ITA18.parameters$h.vec[t]
  R    <- sqrt(dJB^2 + h^2)
  R.mat[t,] <- R
  
  reg.Mlow[t,]  <- ifelse(MAG<=Mh, MAG - Mh, 0)
  reg.Mhigh[t,] <- ifelse(MAG>Mh, MAG - Mh, 0)
  reg.D1[t,]    <- (MAG - Mref)*log10(R)
  reg.D2[t,]    <- log10(R)
  reg.D3[t,]    <- R
}

reg.SS      <- ifelse(SoF=="SS", 1, 0)
reg.TF      <- ifelse(SoF=="TF", 1, 0)
reg.S       <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))

# Create the fd objects by projecting the regressors over an appropriate functional basis

# reg.Mlow
matplot(t.points, reg.Mlow, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.Mlow, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mlow.fd <- smooth.basis(t.points, reg.Mlow, fPar)$fd
plot(Mlow.fd)
dev.off()

# reg.Mhigh
matplot(t.points, reg.Mhigh, type='l')

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mhigh.fd <- smooth.basis(t.points, reg.Mhigh, fPar)$fd
plot(Mhigh.fd)
dev.off()

# R
matplot(t.points, R.mat, type='l', main='True R')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, R.mat, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
R.fd <- smooth.basis(t.points, R.mat, fPar)$fd
plot(R.fd)
title(main='Smoothed R')
dev.off()

# Look what happens for short and long distances: indexes 21 - 36
r <- R.mat[,21]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')
dev.off()

r <- R.mat[,36]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')
dev.off()

# reg.D1
matplot(t.points, reg.D1, type='l')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D1, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D1.fd <- smooth.basis(t.points, reg.D1, fPar)$fd
plot(D1.fd)
dev.off()

# reg.D2
matplot(t.points, reg.D2, type='l')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D2, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D2.fd <- smooth.basis(t.points, reg.D2, fPar)$fd
plot(D2.fd)
dev.off()

# reg.D3
matplot(t.points, reg.D3, type='l')
dev.off()

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D3, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D3.fd <- smooth.basis(t.points, reg.D3, fPar)$fd
plot(D3.fd)
dev.off()

## Construct constant covariates
intercept <- rep(1,n)
xlist <- list(intercept, Mlow.fd, Mhigh.fd, reg.SS, reg.TF, D1.fd, D2.fd, D3.fd, reg.S)
save(xlist, file='output/preprocessed-data/xlist-logg.RData')


## COLLINEARITY ANALYSIS -------------------------------------------------------
#### Load Data -------------------------------------------------------------------
load('output/preprocessed-data/data.RData')
load('output/preprocessed-data/xlist-logg.RData')

t.points <- log10(T.period)
t.points[1] <- -2.5

#### Analysis of multicollinearity in the covariates -----------------------------
x7 <- xlist[[7]]
x8 <- xlist[[8]]
cor.mat.78 <- cor.fd(t.points, x7, t.points, x8)

x6 <- xlist[[6]]
x3 <- xlist[[3]]
cor.mat.36 <- cor.fd(t.points, x3, t.points, x6)

x2 <- xlist[[2]]
cor.mat.23 <- cor.fd(t.points, x2, t.points, x3)

cor.mat.26 <- cor.fd(t.points, x2, t.points, x6)

cor.mat.37 <- cor.fd(t.points, x3, t.points, x7)

xticks <- seq(-2,1,1)
ext.ticks <- c(-2.5,xticks)

levs <- 5

pdf(file = "output/supplementary/geom-and-anel.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.78,
        xlab="Geometric Attenuation",
        ylab="Anelastic Attenuation",
        #color.palette = function(n) hcl.colors(n, "RdYlBu", rev = TRUE),
        main=paste("(a) Correlation across periods for\n",
                   "Geometric and Anelastic Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels = levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
dev.off()

pdf(file = "output/supplementary/highmag-and-geom.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.36,
        xlab="High Magnitudes",
        ylab="Geometric Attenuation",
        main=paste("(d) Correlation across periods for\n",
                   "High Magnitudes and Geometric Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()

pdf(file = "output/supplementary/lowmag-and-geom.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.26,
        xlab="Low Magnitudes",
        ylab="Geometric Attenuation",
        main=paste("(c) Correlation across periods for\n",
                   "Low Magnitudes and Geometric Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()

pdf(file = "output/supplementary/lowmag-and-highmag.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.23,
        xlab="Low Magnitudes",
        ylab="High Magnitudes",
        main=paste("(b) Correlation across periods for\n",
                   "Low and High Magnitudes"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()


