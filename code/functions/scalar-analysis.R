#' Functions to perform a scalar analysis and get the results from Lanzano et al. (2018),
#' ultimately for comparison with the results from the functional analysis.
#' 
#' 


fit_ITA18 <- function(data.f, curves)
{
  # Data -----------------------------------------------------------------------
  dJB      <- data.f$dJB
  MAG      <- data.f$MAG
  SoF      <- data.f$SoF
  VS30     <- data.f$VS30
  
  q <- 9
  
  # Mh, Mref, h ----------------------------------------------------------------
  load('data/ITA18_parameters.RData')
  Mh.vec   <- ITA18.parameters$Mh.vec
  Mref.vec <- ITA18.parameters$Mref.vec
  h.vec    <- ITA18.parameters$h.vec
  
  n <- dim(curves)[2]
  N <- dim(curves)[1]
  
  coefs.ITA18 <- as.data.frame(matrix(data=0, nrow=N, ncol=q))
  y.hat.ITA18 <- as.data.frame(matrix(data=NA, nrow=N, ncol=n))
  sigma.ITA18 <- numeric(N)
  
  # for every period, repeat
  for(t in 1:N) #t=1
  {
    y        <- curves[t,]
    y.in     <- (!is.na(y))
    n.in     <- sum(y.in)
    
    M.h.T    <- Mh.vec[t]
    M.ref.T  <- Mref.vec[t]
    h.T      <- h.vec[t]
    
    # generate dataset for the ITA18 model at time t
    R           <- sqrt(dJB[y.in]^2 + h.T^2)
    reg.M_l     <- ifelse(MAG[y.in]<=M.h.T, MAG[y.in] - M.h.T, 0)
    reg.M_h     <- ifelse(MAG[y.in]>M.h.T, MAG[y.in] - M.h.T, 0)
    reg.SS      <- ifelse(SoF[y.in]=="SS", 1, 0)
    reg.TF      <- ifelse(SoF[y.in]=="TF", 1, 0)
    reg.D_1     <- (MAG[y.in] - M.ref.T)*log10(R)
    reg.D_2     <- log10(R)
    reg.D_3     <- R
    reg.S       <- ifelse(VS30[y.in]<=1500, log10(VS30[y.in]/800), log10(1500/800))
    
    mod.ITA18.T <- lm(y[y.in] ~ reg.M_l + reg.M_h + reg.SS + reg.TF + reg.D_1 + reg.D_2 +
                        reg.D_3 + reg.S)
    
    coefs.ITA18[t,] <- mod.ITA18.T$coefficients
    y.hat.ITA18[t,y.in] <- mod.ITA18.T$fitted.values
    
  }
  names(coefs.ITA18) <- c('a0', 'b1', 'b2', 'f1', 'f2', 'c1', 'c2', 'c3', 'k0')
  
  outlist <- list(coefs.ITA18 = coefs.ITA18,
                  y.hat.ITA18 = y.hat.ITA18)
  
  return(outlist)
}


predict_ITA18 <- function(data, coefs.ITA18)
{
  
  # Data ----------------------------------------------------------------------
  dJB      <- data$dJB
  MAG      <- data$MAG
  SoF      <- data$SoF
  VS30     <- data$VS30
  
  # ITA18 regressors -------------------------------------------------
  a0       <- coefs.ITA18$a0
  b1       <- coefs.ITA18$b1
  b2       <- coefs.ITA18$b2
  c1       <- coefs.ITA18$c1
  c2       <- coefs.ITA18$c2
  c3       <- coefs.ITA18$c3
  k0       <- coefs.ITA18$k0
  f1       <- coefs.ITA18$f1
  f2       <- coefs.ITA18$f2
  
  load('data/ITA18_parameters.RData')
  
  Mh.vec   <- ITA18.parameters$Mh.vec
  Mref.vec <- ITA18.parameters$Mref.vec
  h.vec    <- ITA18.parameters$h.vec
  
  n <- length(dJB)
  N <- length(a0)
  
  SA.ITA18 <- matrix(data=0, nrow=N, ncol=n)
  
  # for every period, repeat
  for(t in 1:N) #t=1
  {
    
    a0.ITA18 <- a0[t]
    b1.ITA18 <- b1[t]
    b2.ITA18 <- b2[t]
    c1.ITA18 <- c1[t]
    c2.ITA18 <- c2[t]
    c3.ITA18 <- c3[t]
    k0.ITA18 <- k0[t]
    f1.ITA18 <- f1[t]
    f2.ITA18 <- f2[t]
    
    beta.T.ITA18 <- as.matrix(c(a0.ITA18, b1.ITA18, b2.ITA18, f1.ITA18,
                                f2.ITA18, c1.ITA18, c2.ITA18, c3.ITA18, k0.ITA18))
    
    M.h.T    <- Mh.vec[t]
    M.ref.T  <- Mref.vec[t]
    h.T      <- h.vec[t]
    
    # generate dataset for the ITA18 model at time t
    R           <- sqrt(dJB^2 + h.T^2)
    intercept   <- rep(1, n)
    reg.M_l     <- ifelse(MAG<=M.h.T, MAG - M.h.T, 0)
    reg.M_h     <- ifelse(MAG>M.h.T, MAG - M.h.T, 0)
    reg.SS      <- ifelse(SoF=="SS", 1, 0)
    reg.TF      <- ifelse(SoF=="TF", 1, 0)
    reg.D_1     <- (MAG - M.ref.T)*log10(R)
    reg.D_2     <- log10(R)
    reg.D_3     <- R
    reg.S       <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))
    
    XX.T        <- cbind(intercept, reg.M_l, reg.M_h, reg.SS, reg.TF, reg.D_1, reg.D_2, reg.D_3, reg.S)
    
    SA.ITA18[t,] <- XX.T%*%beta.T.ITA18
    
  }
  
  return(SA.ITA18)
}