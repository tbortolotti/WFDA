generate_data <- function(seed,
                          case=c("CASE-1", "CASE-2"),
                          reg.info,
                          case.info)
{
  
  ## Utilities
  set.seed(seed)

  t.points       <- reg.info$t.points
  xlist          <- reg.info$xlist
  beta_estimates <- reg.info$beta_estimates
  res            <- reg.info$res
  
  n.sim          <- case.info$n.sim
  if(case=="CASE-1")
  {
    smooth.noise <- case.info$noise
    perc.po      <- case.info$perc
    left.bound   <- case.info$left.bound
    right.bound  <- case.info$right.bound
  } else if (case=="CASE-2"){
    smooth.noise <- case.info$noise
    ext.noise    <- case.info$ext.noise
    perc.noise   <- case.info$perc
    noise.bound  <- case.info$left.bound
    right.bound  <- case.info$right.bound
  }
  
  ## Sample the scalar covariate: x1
  K.fd <- xlist[[2]]
  K.val <- eval.fd(1, K.fd)
  mu1 <- mean(K.val)
  sd1 <- sd(K.val)
  #x1 <- rnorm(n.sim, mu1, sd1)
  #x1 <- rnorm(n.sim, 10*mu1, sd1)
  x1 <- rnorm(n.sim, 0, 2.5*sd1)
  
  ## Sample the functional covariate: x2
  x2.fd <- xlist[[3]]
  x2.fd$coefs <- xlist[[3]]$coefs[,1:n.sim]
  
  ## expand the coefficients of a random term between 1 and 1.2
  fact.to.exp <- matrix(data=runif(dim(xlist[[3]]$coefs)[1]*n.sim, min=1, max=1.2),
                        ncol=n.sim)
  expanded.coefs <- xlist[[3]]$coefs[,1:n.sim]*fact.to.exp
  mu.vec <- as.numeric(apply(X=expanded.coefs, MARGIN=1, FUN=mean))
  cov.mat <- var(t(expanded.coefs))
  x2.fd$coefs <- t(mvrnorm(n.sim, mu.vec, cov.mat))
  
  ## FPCA of the residuals
  pca <- pca.fd(res, nharm=6, centerfns=TRUE)
  
  harm1 <- pca$harmonics[1]
  harm1$coefs <- harm1$coefs %*% rep(1,n.sim)
  
  harm2 <- pca$harmonics[2]
  harm2$coefs <- harm2$coefs %*% rep(1,n.sim)
  
  pca.scores <- pca$scores[,1:2]
  scores.mu.vec <- as.numeric(apply(X=pca.scores, MARGIN=2, FUN=mean))
  scores.cov.mat <- var(pca.scores)
  new.scores <- mvrnorm(n.sim, scores.mu.vec, scores.cov.mat)
  
  res_tilde <- new.scores[,1]*harm1 + new.scores[,2]*harm2
  
  # xlist for WFDA
  x1 <- as.matrix(x1)
  intercept <- as.matrix(rep(1,n.sim))
  xlist.new <- list(intercept, x1, x2.fd)
  
  ## Construction of y
  beta0 <- beta_estimates[[1]]$fd
  beta0$coefs <- beta0$coefs/3
  beta1 <- beta_estimates[[2]]$fd
  beta2 <- beta_estimates[[3]]$fd
  
  for(i in 1:n.sim) #i=1
  {
    x2i <- x2.fd
    x2i$coefs <- as.matrix(x2.fd$coefs[,i], ncol=1)
    
    resi <- res_tilde
    resi$coefs <- as.matrix(res_tilde$coefs[,i], ncol=1)
    
    yi <- beta0 + x1[i]*beta1 + x2i*beta2 + resi
    yi.nores <- beta0 + x1[i]*beta1 + x2i*beta2
    # yi.nointerc <- x1[i]*beta1 + x2i*beta2
    
    if(i==1)
    {
      y.fd <- yi
      y.fd$coefs <- yi$coefs %*% rep(1,n.sim)
      
      y.fd.nores <- yi.nores
      y.fd.nores$coefs <- yi.nores$coefs %*% rep(1,n.sim)
      
      # y.fd.nointerc <- yi.nointerc
      # y.fd.nointerc$coefs <- yi.nointerc$coefs %*% rep(1,n.sim)
    }
    
    y.fd$coefs[,i] <- yi$coefs
    y.fd.nores$coefs[,i] <- yi.nores$coefs
    # y.fd.nointerc$coefs[,i] <- yi.nointerc$coefs
    
  }
  # 
  # par(mfrow=c(2,3))
  # plot(beta1, ylab='value')
  # plot(beta2, ylab='value')
  # plot(beta0, ylab='value')
  # plot(x1, ylab='values')
  # plot(x2.fd, ylab='values')
  # plot(y.fd, ylab='values')
  # 
  # for(i in 1:n.sim) #i=1
  # {
  #   x2i <- x2.fd
  #   x2i$coefs <- as.matrix(x2.fd$coefs[,i], ncol=1)
  # 
  #   x1effect <- x1[i]*beta1
  #   x2effect <- x2i*beta2
  #   if(i==1)
  #   {
  #     x1effect.fd <- x1effect
  #     x1effect.fd$coefs <- x1effect$coefs %*% rep(1,n.sim)
  # 
  #     x2effect.fd <- x2effect
  #     x2effect.fd$coefs <- x2effect$coefs %*% rep(1,n.sim)
  #   }
  # 
  #   x1effect.fd$coefs[,i] <- x1effect$coefs
  #   x2effect.fd$coefs[,i] <- x2effect$coefs
  # }
  # 
  # par(mfrow=c(2,2))
  # plot(beta0, ylab='value')
  # plot(x1effect.fd, ylab='value')
  # plot(x2effect.fd, ylab='value')
  # plot(y.fd.nores, ylab='values')
  # 
  # plot(y.fd.nointerc, ylab='values')
  
  ## Add the smoothing errors
  ## Add an iid error to the curves
  t.grid <- seq(range(t.points)[1], range(t.points)[2], 0.25)
  
  if(case == "CASE-1")
  {
    ## CASE 1
    ## Partially observed functional data
    Ui <- runif(n.sim, min=left.bound,max=right.bound)
    Pi <- rbinom(n.sim, size=1, prob=perc.po)
    Pi.logic <- ifelse(Pi==1, TRUE, FALSE)
    T_hp <- rep(right.bound, n.sim)
    T_hp[Pi.logic] <- Ui[Pi.logic]
    
    y.noisy <- eval.fd(t.grid, y.fd)
    for(i in 1:n.sim)
    {
      y.noisy[,i] <- y.noisy[,i] + rnorm(length(t.grid), 0, smooth.noise)
    }
    
    y.no.noise <- eval.fd(t.grid, y.fd)
    y.po <- y.noisy
    for(i in 1:n.sim)
    {
      t.idx <- (t.grid>T_hp[i])
      y.po[t.idx,i] <- NA
    }
    
    
    ## OUTPUT
    out <- list(t.points        = t.grid,
                T_hp            = T_hp,
                curves          = y.po,
                curves.true     = y.noisy,
                curves.true.fd  = y.fd.nores,
                curves.reg.fd   = y.fd,
                xlist           = xlist.new)
    return(out)
  }
  
}