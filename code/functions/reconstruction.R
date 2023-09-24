
#' Reconstruction method considered in this work:
#' - Kraus (2015)
#' - Kneip and Liebl ()
#' - extrapolation
#' 
#' The functions for the methods by Kraus (2015) and by Kneip and Liebl (2020)
#' are those given in the reconstPoFD package of Kneip and Liebl (2020).
#' As we need the estimated variability of the reconstructed trajectory to define
#' reconstruction-driven functional weights, the original functions are reported here
#' with a minor modification that allows one to estimate such variability and return
#' it as part of the output.

## KRAUS ----------------------------------------------------------------------


#' Reconstruct functional data as proposed by Kraus (JRSSB, 2015)
#'
#' This function allows to reconstruct partially observed functional data as proposed in: 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' @param X_mat        pxn matrix (p: number of discretization points, n=number of functions
#' @param alpha        Ridge parameter. If alpha = NULL (default), an optimal alpha is determined by GCV
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed. Default (reconst_fcts=NULL) will reconstruct all functions.
#' @export reconstructKraus
#' @examples  
#' SimDat       <- simuldata(n = 50, a = 0, b = 1, DGP="DGP3")
#' Y_mat        <- SimDat[['Y_mat']]
#' U_mat        <- SimDat[['U_mat']]
#' U_true_mat   <- SimDat[['U_true_mat']]
#' ##
#' result        <- reconstructKraus(X_mat = Y_mat)
#' Y_reconst_mat <- result[['X_reconst_mat']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", 
#' main="Original Data", ylab="", xlab="")
#' matplot(x=U_true_mat[,1:5], y=Y_reconst_mat[,1:5], col=gray(.5), 
#' type="l", main="Kraus (2015)", ylab="", xlab="")
#' par(mfrow=c(1,1))
reconstructKraus <- function(X_mat, 
                             alpha        = NULL,
                             reconst_fcts = NULL){
  ##
  mean_vec      <- meanKraus(X_mat)
  cov_mat       <- covKraus(X_mat)
  n             <- ncol(X_mat)
  if(is.null(reconst_fcts)){
    reconst_fcts <- 1:n
  }
  X_reconst_mat <- X_mat[, reconst_fcts, drop=FALSE]
  W_reconst_mat <- matrix(data=1, nrow=dim(X_reconst_mat)[1], ncol=length(reconst_fcts))
  ##
  NonNA_fcts    <- apply(X_mat,2,function(x)!any(is.na(x)))
  X_Compl_mat   <- X_mat[,NonNA_fcts]
  alpha_vec     <- rep(NA, length(reconst_fcts)) 
  df_vec        <- rep(NA, length(reconst_fcts)) 
  ##
  for(i in 1:length(reconst_fcts)){ #i=1
    X_tmp      <- X_mat[,reconst_fcts[i]]
    ##
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    ##
    if(is.null(alpha)){
      alpha_vec[i] <- stats::optimize(f = function(alpha){gcvKraus(cov_mat     = cov_mat, 
                                                                   mean_vec    = mean_vec, 
                                                                   X_Compl_mat = X_Compl_mat, 
                                                                   M_bool_vec  = M_bool_vec, 
                                                                   alpha       = alpha)},
                                      interval = c(.Machine$double.eps, sum(diag(cov_mat))*n), maximum = FALSE)$minimum
    }else{
      alpha_vec[i] <- alpha
    }
    ##
    result_tmp <- reconstKraus_fun(cov_mat    = cov_mat, 
                                   X_cent_vec = c(X_tmp - mean_vec), 
                                   alpha      = alpha_vec[i])
    ##
    X_reconst_mat[,i]  <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    df_vec[i]          <- result_tmp[['df']]
    W_reconst_mat[M_bool_vec,i]  <- 1-result_tmp[['hi']]
  }
  return(list("X_reconst_mat" = X_reconst_mat,
              "alpha"         = alpha_vec, 
              "df"            = df_vec,
              "W_reconst_mat" = W_reconst_mat))
}

meanKraus <- function(X_mat){
  rowMeans(X_mat, na.rm = TRUE)
}

covKraus <- function(X_mat){
  p <- nrow(X_mat)
  n <- ncol(X_mat)
  ##
  X_cent_mat  <- X_mat - rowMeans(X_mat, na.rm = TRUE)
  ##
  covKraus_mat <- matrix(NA, ncol = p, nrow = p)	
  for(s in seq(1, p)){
    for(t in seq(s, p)){
      X_cent_s  <- X_cent_mat[s, ]
      X_cent_t  <- X_cent_mat[t, ]
      n_na      <- sum(is.na(c(X_cent_s * X_cent_t)))
      if(n-n_na == 0){
        covKraus_mat[s,t] <- NA
      }else{
        covKraus_mat[s,t] <- mean(X_cent_s * X_cent_t, na.rm = TRUE)
      }
      covKraus_mat[t,s] <- covKraus_mat[s,t]
    }
  }
  return(covKraus_mat)
}

reconstKraus_fun <- function(cov_mat, X_cent_vec, alpha=1e-4){
  ##
  M_bool_vec       <- is.na(X_cent_vec)
  O_bool_vec       <- !M_bool_vec
  p                <- nrow(cov_mat)
  ##
  covMO_mat        <- cov_mat[M_bool_vec, O_bool_vec] 
  covOO_mat        <- cov_mat[O_bool_vec, O_bool_vec]
  ##
  covOO_a_mat      <- covOO_mat + alpha * diag(p)[O_bool_vec, O_bool_vec]
  covOO_a_mat_inv  <- solve(covOO_a_mat)
  Aa_mat           <- covMO_mat %*% covOO_a_mat_inv
  ##
  X_M_fit_cent_vec <- Aa_mat %*% X_cent_vec[O_bool_vec]
  ##
  X_cent_reconst_vec             <- X_cent_vec
  X_cent_reconst_vec[M_bool_vec] <- X_M_fit_cent_vec
  ## compute df
  #df <- sum(diag(covOO_a_mat_inv %*% covOO_mat))
  lam00 <- base::eigen(covOO_mat)$values
  lam00 <- lam00[lam00>0]
  df    <- sum(lam00/(lam00+alpha))
  ##
  
  ## compute h(t)
  covMM_mat <- cov_mat[M_bool_vec, M_bool_vec]
  V_mat <- covMM_mat - Aa_mat %*% covOO_mat %*% t(Aa_mat)
  if(all(diag(V_mat)>0)){
    v2i <- diag(V_mat)
    vi <- sqrt(v2i)
    h0 <- 0.2 * max(vi)
    hi <- ifelse(h0<vi, vi, h0)
  }else{
    hi <- rep(0.01,sum(M_bool_vec))
  }
  if(length(hi)==1){
    hi_scaled <- hi/sqrt(covMM_mat)
  }else{
    hi_scaled <- hi/sqrt(diag(covMM_mat))
  }
  
  return(list("X_cent_reconst_vec" = X_cent_reconst_vec,
              "df"                 = df,
              "alpha"              = alpha,
              "hi"                 = hi_scaled))
}

gcvKraus <- function(cov_mat, mean_vec, X_Compl_mat, M_bool_vec, alpha){
  n_Compl  <- ncol(X_Compl_mat)
  rss_vec  <- rep(NA,n_Compl)
  ##
  for(j in 1:n_Compl){
    X_gcv             <-  X_Compl_mat[,j]
    X_gcv[M_bool_vec] <- NA
    ##
    result_tmp <- reconstKraus_fun(cov_mat    = cov_mat, 
                                   X_cent_vec = c(X_gcv - mean_vec), 
                                   alpha      = alpha)
    ##
    X_fit      <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    ##
    rss_vec[j] <- sum((X_fit[M_bool_vec] - X_Compl_mat[M_bool_vec,j])^2)
  }
  gcv <- sum(rss_vec)/((1-result_tmp[['df']]/n_Compl)^2)
  ##
  return(gcv)
}
gcvKraus <- Vectorize(FUN = gcvKraus, vectorize.args = "alpha")


## KNEIP AND LIEBL -----------------------------------------------------------

#' Reconstruct partially observed functions
#'
#' This function allows you to reconstruct the missing parts of a function given the observed parts.
#' @param Ly           List of Y-values. The ith (i=1,...,n) list-element contains \eqn{Y_{i1},\dots,Y_{im}}{Y_{i1},...,Y_{im}}
#' @param Lu           List of U-values. The ith (i=1,...,n) list-element contains \eqn{U_{i1},\dots,U_{im}}{U_{i1},...,U_{im}}
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed. Default (reconst_fcts=NULL) will reconstruct all functions.
#' @param method       One of the following options: 'Error=0_AlignYES_CommonGrid', 'Error>0_AlignYES', 'Error>=0_AlignNO', 'Error>0_AlignYES_CEscores', 'Error>0_AlignNO_CEscores', 'PACE', 'Error=0_PACE'.
#' @param K            Truncation parameter. If K=NULL (default), K is determined using GCV.
#' @param nRegGrid     Number of gridpoints within [a,b] used for the reconstruction result.
#' @param maxbins      If maxbins=NULL (default), maxbins is set to 1000. For speeding up simulations, use, for instance maxbins=200.
#' @param progrbar     Show progress bar (TRUE) or not (FALSE, default)
#' @export reconstructKneipLiebl
#' @examples  
#' 
#' a <- 0; b <- 1
#' set.seed(223109)
#' 
#' ## Generate partially observed functional data with error
#' SimDat        <- simuldata(n = 50, m=15, a = a, b = b, DGP="DGP1")
#' ## 
#' Y_list   <- SimDat[['Y_list']]; Y_mat <- SimDat[['Y_mat']]
#' U_list   <- SimDat[['U_list']]; U_mat <- SimDat[['U_mat']]
#' ##
#' ## Reconstruction with alignments of reconstructed parts and presmoothed fragments
#' ## using integral scores
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignYES', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction with out alignments
#' ## using integral scores
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>=0_AlignNO', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))
#' 
#' 
#' ## Reconstruction with alignments of reconstructed parts and presmoothed fragments
#' ## using conditional expectation scores
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignYES_CEscores', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction with out alignments
#' ## using conditional expectation scores
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignNO_CEscores', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))
#' 
#' 
#' ## Generate partially observed functional data without error 
#' SimDat        <- simuldata(n = 50, a = a, b = b, DGP="DGP2")
#' ## 
#' Y_list   <- SimDat[['Y_list']]; Y_mat <- SimDat[['Y_mat']]
#' U_list   <- SimDat[['U_list']]; U_mat <- SimDat[['U_mat']]
#' ##
#' ## Reconstruction with alignments of reconstructed parts and observed fragments
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error=0_AlignYES_CommonGrid', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction without alignments
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>=0_AlignNO', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))

reconstructKneipLiebl <- function(Ly,
                                  Lu, 
                                  reconst_fcts = NULL,
                                  method       = c('Error=0_AlignYES_CommonGrid',
                                                   'Error>0_AlignYES',
                                                   'Error>=0_AlignNO',
                                                   'Error>0_AlignYES_CEscores',
                                                   'Error>0_AlignNO_CEscores',
                                                   'PACE',
                                                   'Error=0_PACE'),
                                  K            = NULL,
                                  nRegGrid     = NULL,
                                  maxbins      = NULL, 
                                  progrbar     = FALSE){
  ##
  method <- switch(method, 
                   `Error=0_AlignYES_CommonGrid` = 1, 
                   `Error>0_AlignYES`            = 2,
                   `Error>=0_AlignNO`            = 3,
                   `Error>0_AlignYES_CEscores`   = 4,
                   `Error>0_AlignNO_CEscores`    = 5,
                   PACE                          = 6,
                   `Error=0_PACE`                = 7)
  ##
  center <- TRUE
  ##
  n      <- length(Ly)
  if(is.null(reconst_fcts)){ reconst_fcts <- 1:n }
  ##
  Y_reconst_list <- vector("list", length(reconst_fcts))
  U_reconst_list <- vector("list", length(reconst_fcts))
  W_reconst_list <- vector("list", length(reconst_fcts))
  K_vec          <- rep(NA,        length(reconst_fcts))
  ##
  if(method == 1){# 'Error=0_AlignYES_CommonGrid'
    ## method 1: 
    ## Requires: error=0, common argvalues (i.e., data structure as in Kraus JRSSB)
    ## Uses:     Classical (integral) scores. Alignment of the reconstructed parts and the fully observed fragments.
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( any(fpca_obj$obs_argvalsO[[i]] != fpca_obj$argvalsO[[i]]) ){stop("The fragment must be fully observed (obs_argvalsO == argvalsO).") }
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 1,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      fragmO     <- c(stats::na.omit(c(fpca_obj$Y[i,])))
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          evaluesO    = fpca_obj$evaluesO[[i]],
                                          fragmO      = fragmO, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 2){# 'Error>0_AlignYES'
    ## method 2: 
    ## Requires: error>0, random or common argvalues (needs a large number of observed argvalues per function)
    ## Does:     Classical (integral) scores. Pre-smoothing of observed fragments. Alignment of reconstructed parts and pre-smoothed fragments.
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( !all(range(fpca_obj$obs_argvalsO[[i]]) == range(fpca_obj$argvalsO[[i]])) ){
        stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO.") }
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 2,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit(c(fpca_obj$Y[i,]))), x=fpca_obj$obs_argvalsO[[i]]))
      fragmO_presmooth  <- stats::predict(smooth.fit, fpca_obj$argvalsO[[i]])$y
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          fragmO      = fragmO_presmooth, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 3){# 'Error>=0_AlignNO'
    ## method 3: 
    ## Requires: error>=0, random or common argvalues
    ## Uses:     Classical (integral) scores. No alignment 
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 3,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          fragmO      = NULL, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 4){# 'Error>0_AlignYES_CEscores'
    ## method 4: 
    ## Requires: error>0, random or common argvalues (needs a large number of observed argvalues per function for the pre-smoothing)
    ## Uses:     CEscores (PACE). Pre-smoothing of observed fragments. Alignment of reconstructed parts with the pre-smoothed fragments.
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( !all(range(fpca_obj$obs_argvalsO[[i]]) == range(fpca_obj$argvalsO[[i]])) ){
        stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO.") }
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 4,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit(c(fpca_obj$Y[reconst_fcts[i],]))), x=fpca_obj$obs_argvalsO[[i]]))
      fragmO_presmooth  <- stats::predict(smooth.fit, fpca_obj$argvalsO[[i]])$y
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu,
                                          cov         = fpca_obj$cov,
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$CE_scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          evaluesO    = fpca_obj$evaluesO[[i]],
                                          fragmO      = fragmO_presmooth, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
      W_reconst_list[[i]]   <- result_tmp[['w_reconst']]
    }
  }
  if(method == 5){# 'Error>0_AlignNO_CEscores'
    ## method 5: 
    ## Requires: error>0, random or common argvalues (works also for a few observed argvalues per function, since no pre-smoothing)
    ## Uses:     CEscores (PACE). No alignment.
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 2
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 5,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          cov         = fpca_obj$cov,
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$CE_scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          evaluesO    = fpca_obj$evaluesO[[i]],
                                          fragmO      = NULL, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
      W_reconst_list[[i]]  <- result_tmp[['w_reconst']]
    }
  }
  if(method == 6){# 'PACE'
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        PACE         = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 6,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- PACE_fun(mu          = fpca_obj$mu, 
                             argvals     = fpca_obj$argvals, 
                             scoresP     = fpca_obj$scoresP[[i]], 
                             efunctionsP = fpca_obj$efunctionsP, 
                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
    }
  } 
  if(method == 7){# 'PACE_Error==0'
    fpca_obj <- my.fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        PACE_E0      = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 7,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- PACE_fun(mu          = fpca_obj$mu, 
                             argvals     = fpca_obj$argvals, 
                             scoresP     = fpca_obj$scoresP[[i]], 
                             efunctionsP = fpca_obj$efunctionsP, 
                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
    }
  }
  ##
  if(!is.null(nRegGrid)){
    ## Evaluate the reconstruced functions at a regular gird of length nRegGrid
    xout <- seq(from = min(fpca_obj$argvals), to = max(fpca_obj$argvals), len=nRegGrid)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if(all(is.na(Y_reconst_list[[i]]))){
        Y_reconst_list[[i]]  <- rep(NA, length(xout))
        U_reconst_list[[i]]  <- xout
      } else {
        Reconstr_on_RegGrid  <- stats::spline(y = Y_reconst_list[[i]], x = U_reconst_list[[i]], xout = xout)
        Y_reconst_list[[i]]  <- Reconstr_on_RegGrid$y
        U_reconst_list[[i]]  <- Reconstr_on_RegGrid$x
      }
    }
  }
  ##
  return(list(
    "Y_reconst_list"  = Y_reconst_list,
    "U_reconst_list"  = U_reconst_list,
    "W_reconst_list"  = W_reconst_list,
    "K"               = K_vec
  ))
}


gcvKneipLiebl <- function(fpca_obj, argvalsO, method, pev = 0.99, progrbar = FALSE){
  ##
  Y           <- fpca_obj$Y # data pre-processed by irreg2mat()
  mu          <- fpca_obj$mu
  cov         <- fpca_obj$cov
  sigma2      <- fpca_obj$sigma2
  argvals     <- fpca_obj$argvals
  efunctionsP <- fpca_obj$efunctionsP
  evaluesP    <- fpca_obj$evaluesP
  npcP        <- length(evaluesP)
  a           <- min(argvals)
  b           <- max(argvals)
  ##
  ## a function is considered as complete if: 
  ## 1. its first and its last observation are non-missing
  ## 2. not considered so far
  compl_fcts <- apply(Y, 1, function(x){
    !is.na(utils::head(x, 1)) & !is.na(utils::tail(x, 1)) 
    # & length(c(stats::na.omit(x)))>=floor(.8*length(argvals))
  })
  ##
  # functions to be reconstructed
  Y.compl       <- Y[compl_fcts,,drop=FALSE]
  # make artifical fragements (all over argvalsO)
  locO          <- match(argvalsO, argvals)
  locM          <- c(1:length(argvals))[-locO]
  Y.pred        <- Y.compl
  Y.pred[,locM] <- NA
  ## the artifical fragments must have sufficiently many observations (>=5) for the pre-smoothing
  slct          <- apply(Y.pred, 1, function(x){length(c(stats::na.omit(x)))>=5})
  Y.pred        <- Y.pred[slct,,drop=FALSE]
  n_compl       <- nrow(Y.pred)
  ##
  if(any(method == c(1,2,3,4,5))){ # all, but PACE and PACE_E0
    muO               <- mu[locO]
    w                 <- quadWeights(argvalsO, method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    # CovOO
    V                 <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalues           <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues           <- replace(evalues, which(evalues <= 0), 0)
    ##
    npc               <- length(evalues[evalues>0])
    npc               <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
    ##
    efunctionsO       <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evaluesO          <- evalues[1:npc]  
    D.inv             <- diag(1/evaluesO, nrow = npc, ncol = npc)
    Z                 <- efunctionsO
    #'# ##################################################################
    ## Reconstructive eigenfunctions
    efun_reconst      <- matrix(NA, nrow=length(argvals), ncol=npc)
    ##
    for(k in seq_len(npc)){
      efun_reconst[,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                FUN = function(x){pracma::trapz(x=argvalsO, efunctionsO[,k] * x)})
      efun_reconst[,k] <- efun_reconst[,k] / evaluesO[k]
    }
    #'# ##################################################################
  } else {
    npc <- npcP  
  }
  ##
  if(n_compl <= 1){
    # warning("Too few complete functions; we use the pve=0.9 criterium")
    if(any(method == c(6,7))){
      K.pve <- which(cumsum(evaluesP[evaluesP>0])/sum(evaluesP[evaluesP>0])>=.9)[1]
    } else {
      K.pve <- which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=.9)[1]
    }
    return(K.pve)
  }
  ##
  rss_mat           <- matrix(NA, n_compl, npc)
  if(progrbar){
    ## Progress Bar
    cat("Select K via GCV:\n")
    pb                <- utils::txtProgressBar(min = 0, max = n_compl*npc, style = 3)
    counter           <- 0
  }
  ##
  for(i in seq_len(n_compl)){# i <- 3
    #'# Needed for all following ###########################################
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, ncol(Y)))
    obs_locO          <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvalsO))
    #'# ####################################################################
    ## 
    if(any(method==c(1,2,3))){
      ## Classical scores (via intergration)
      scoresO      <- apply(X      = efunctionsO[obs_locO,,drop=FALSE], 
                            MARGIN = 2, 
                            FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=argvalsO[obs_locO])})
    }
    ##
    if(any(method==c(2,4))){
      ## Pre-smoothing of fragmO
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit((Y.pred[i,]))), x=argvalsO[obs_locO]))
      fragmO_presmooth  <- stats::predict(smooth.fit, argvalsO)$y
    }
    ##
    if(any(method==c(4,5))){
      ## CEscores (PACE scores with respect to efunctionsO)
      if(sigma2           ==  0){sigma2 <- 1e-6}
      if(length(obs_locO) < npc){
        npcO        <- length(obs_locO)
        Zcur        <- Z[obs_locO,1:npcO,drop=FALSE]
        ZtZ_sD.inv  <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npcO,1:npcO])
        CE_scoresO  <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
        CE_scoresO  <- c(CE_scoresO, rep(0, npc-npcO))
      }else{
        Zcur        <- Z[obs_locO,,drop=FALSE]
        ZtZ_sD.inv  <- solve(crossprod(Zcur) + sigma2 * D.inv)
        CE_scoresO  <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
      }
    }
    if(any(method == c(6,7))){
      ## PACE (PACE scores with respect to efunctionsP)
      obs_locP        <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvals))
      Y.centP         <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, ncol(Y)))
      D.invP          <- diag(1/evaluesP, nrow = npc, ncol = npc)
      ##
      if(sigma2 == 0){sigma2 <- 1e-6}
      if(method == 7){sigma2 <- 0}
      if(length(obs_locP) < npc){
        npcO            <- length(obs_locP)
        ZcurP           <- efunctionsP[obs_locP, 1:npcO, drop=FALSE]
        ##
        ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcO, 1:npcO]), silent = TRUE)
        # while(is.error(ZtZ_sD.invP)){# for preventing singularities
        #   sigma2          <- sigma2 + .Machine$double.eps
        #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcO, 1:npcO]), silent = TRUE)
        # }
        if(is.error(ZtZ_sD.invP)){
          ZtZ_sD.invP     <- matrix(NA, npcO, npcO)
        }
        scoresP         <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
        scoresP         <- c(scoresP, rep(0, npc-npcO))
      }else{
        ZcurP           <- efunctionsP[obs_locP, , drop=FALSE]
        ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npc, 1:npc]), silent = TRUE)
        # while(is.error(ZtZ_sD.invP)){# for preventing singularities
        #   sigma2          <- sigma2 + .Machine$double.eps
        #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npc, 1:npc]), silent = TRUE)
        # }
        if(is.error(ZtZ_sD.invP)){
          ZtZ_sD.invP     <- matrix(NA, npc, npc)
        }
        scoresP         <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
      }
    }
    # if(method == 8){
    #   if(sigma2 == 0){sigma2 <- 1e-6}
    #   ##
    #   efun_reconst_orth <- pracma::gramSchmidt(efun_reconst)$Q
    #   Zcur              <- efun_reconst_orth
    #   Zcur              <- Zcur[obs_locO,1:npc,drop=FALSE]
    #   ZtZ_sD.inv        <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npc,1:npc])
    #   CE_scores_orth    <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    # }
    ##
    for(k in seq_len(npc)){# k <- 5
      ##
      if(method == 1){# Classical scores, with alignment of reconstr and fully observed fragmO
        result_tmp <- reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=Y.pred[i,locO], K=k)
      }
      if(method == 2){# Classical scores, with alignment of reconstr and fragmO
        result_tmp <- reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=fragmO_presmooth, K=k)
      }
      if(method == 3){# Classical scores, without alignment of reconstr and pre-smoothed fragmO
        result_tmp <- reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=NULL,K=k)
      }
      if(method == 4){# CEscores, with alignment of reconstr and pre-smoothed fragmO
        result_tmp <- reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=CE_scoresO,efun_reconst=efun_reconst,fragmO=fragmO_presmooth,K=k)
      }
      if(method == 5){# CEscores, without alignment of reconstr and fragmO
        result_tmp <- reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=CE_scoresO,efun_reconst=efun_reconst,fragmO=NULL,  K=k)
      }
      if(any(method == c(6,7))){# PACE
        result_tmp <- PACE_fun(mu=mu, argvals=argvals, scoresP=scoresP, efunctionsP=efunctionsP, K=k)
      }
      # if(method == 8){# Orth
      #   result_tmp <- reconstKneipLiebl_orth_fun(mu=mu,argvals=argvals,CE_scores_orth=CE_scores_orth,efun_reconst_orth=efun_reconst_orth,K=k)
      # }
      ##
      rss_mat[i,k] <- sum((result_tmp[['y_reconst']][locM] - Y.compl[i,locM])^2, na.rm = TRUE)
      # cat("i=",i," k=",k,"\n")
      if(progrbar){
        counter <- counter+1
        utils::setTxtProgressBar(pb, counter)
      }
    }
  }
  if(progrbar){close(pb)}
  ##
  gcv_k_vec <- colSums(rss_mat)/((1-1:npc/n_compl)^2)
  K.gcv     <- which.min(gcv_k_vec)
  ##
  return(K.gcv)
}

reconstKneipLiebl_fun <- function(mu, cov=NULL, argvals, argvalsO, scoresO, efun_reconst, evaluesO=NULL, fragmO=NULL, K=NULL){
  ##
  K        <- min(ncol(efun_reconst), K)
  locO     <- match(argvalsO, argvals)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(scoresO[1:K, drop=FALSE]) %*% t(efun_reconst[,1:K, drop=FALSE]) ))
  ##
  if( !is.null(fragmO) ){
    ## Aligning the reconstructed parts to fragm0
    reconstr[1:min(locO)]               <- reconstr[1:min(locO)]               + fragmO[1]              - reconstr[min(locO)]
    reconstr[max(locO):length(argvals)] <- reconstr[max(locO):length(argvals)] + fragmO[length(fragmO)] - reconstr[max(locO)]
    reconstr[locO]                      <- fragmO
  }
  ##
  if(is.null(evaluesO) | is.null(cov)){
    w_reconstr <- NULL
  }else{
    idxsM <- (max(locO)+1):length(argvals)
    v2_reconstr <- diag(cov)[idxsM] - (t(evaluesO[1:K, drop=FALSE]) %*% t(efun_reconst[,1:K, drop=FALSE]^2))[idxsM]
    # w_reconstr <- rep(1,length(argvals))
    # w_reconstr[idxsM] <- 1 - sqrt(v_reconstr/diag(cov)[idxsM])
    # if(any(v2_reconstr>0)){
    #   v_reconstr <- sqrt(v2_reconstr[v2_reconstr>0])
    #   h0 <- 0.2*max(v_reconstr)
    #   v_hat_reconstr <- rep(h0,length(idxsM))
    #   v_hat_reconstr[v2_reconstr>0] <- ifelse(h0<v_reconstr,v_reconstr, h0)
    # }else{
    #   v_hat_reconstr <- rep(0.1, length(idxsM))
    # }
    if(all(v2_reconstr>0)){
      v_reconstr <- sqrt(v2_reconstr)
      h0 <- 0.2*max(v_reconstr)
      v_hat_reconstr <- ifelse(h0<v_reconstr,v_reconstr, h0)
    }else{
      v_hat_reconstr <- rep(0.01, length(idxsM))
    }
    w_reconstr <- rep(1,length(argvals))
    w_reconstr[idxsM] <- 1 - v_hat_reconstr/sqrt(diag(cov)[idxsM])
  }
  
  #'# ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals),
              "w_reconst"  = c(w_reconstr)))
  #'# ######################
}

reconstKneipLiebl_orth_fun <- function(mu, argvals, CE_scores_orth , efun_reconst_orth, K=NULL){
  ##
  K        <- min(ncol(efun_reconst_orth), K)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(CE_scores_orth[1:K]) %*% t(efun_reconst_orth[,1:K]) ))
  ##
  #'# ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  #'# ######################
}


PACE_fun <- function(mu, argvals, scoresP, efunctionsP, K=NULL){
  ##
  K        <- min(ncol(efunctionsP), K)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(scoresP[1:K]) %*% t(efunctionsP[,1:K]) ))
  ##
  #'# ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  #'# ######################
}

my.fpca <- function(Ly, Lu, reconst_fcts = NULL, pev = 0.99, CEscores = FALSE, PACE = FALSE, PACE_E0 = FALSE, center = TRUE, maxbins = NULL){
  
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ##
  ydata  <-  data.frame(".id"    = id_vec, 
                        ".index" = unname(unlist(Lu)), 
                        ".value" = unname(unlist(Ly)))
  ##
  nbasis         <- 10
  maxbins        <- ifelse(is.null(maxbins), 1000, maxbins)
  useSymm        <- FALSE # if true, there will be no smoothing accross the diagonal
  makePD         <- FALSE # if true, the NP-estimate of cov is forced to be positive definite
  ##
  # Y (nobs x length(argvals)) with NA's if no data at a argvalue
  Y        <- irreg2mat(ydata, binning = TRUE, maxbins = maxbins) 
  argvals  <- as.numeric(colnames(Y))
  ##
  # functions to be reconstructed
  if(is.null(reconst_fcts)){reconst_fcts <- 1:n}
  Y.pred   <- Y[reconst_fcts,,drop=FALSE]
  # argvals of observed fragments to be reconstructed
  argvalsO <- vector("list", length = length(reconst_fcts))
  for(i in seq_len(length(reconst_fcts))){
    minmax        <- range(as.numeric(names(c(stats::na.omit(Y.pred[i,])))))
    argvalsO[[i]] <- argvals[argvals >= minmax[1] & argvals <= minmax[2]]
  }
  ##
  D      <- NCOL(Y)      # nobs per function
  I      <- NROW(Y)      # number of functions
  I.pred <- NROW(Y.pred) # number of functions to be reconstruced
  ##
  d.vec <- rep(argvals, each = I)
  id    <- rep(1:I, rep(D, I))
  
  #'# MEAN ##############################################################
  if(center){
    ## mean
    gam0    <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu      <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
  }else{
    Y.tilde <- Y
    mu      <- rep(0, D)
  }
  #'# plot(x=argvals,y=mu, type="l")
  #'#
  
  #'# COV ###############################################################
  ## 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
  ## 2. smooth this matrix
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
    cov.sum[  obs.points, obs.points] <- cov.sum[  obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0   <- diag(G.0)
  diag(G.0) <- NA
  if(!useSymm){
    row.vec <- rep(argvals, each = D)
    col.vec <- rep(argvals, D)
    npc.0   <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                                  weights = as.vector(cov.count)), 
                                        newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
    npc.0 = (npc.0 + t(npc.0))/2
  }else{
    use          <- upper.tri(G.0, diag = TRUE)
    use[2, 1]    <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
    usecov.count <- cov.count
    usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
    usecov.count <- as.vector(usecov.count)[use]
    use          <- as.vector(use)
    vG.0         <- as.vector(G.0)[use]
    row.vec      <- rep(argvals, each = D)[use]
    col.vec      <- rep(argvals, times = D)[use]
    mCov         <- mgcv::gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
    npc.0        <- matrix(NA, D, D)
    spred        <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
    tpred        <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
    # Estimated covariance function:
    smVCov       <- mgcv::predict.gam(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
    npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
    npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    # slct <- seq.int(1,length(argvals),len=25)
    # persp(z=npc.0[slct,slct],x=argvals[slct],y=argvals[slct])
  }
  if(makePD){
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
      as.matrix(tmp$mat)
    }
  }
  # Nnumerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  w          <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt      <- diag(sqrt(w))
  Winvsqrt   <- diag(1/(sqrt(w)))
  V          <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues    <- replace(evalues, which(evalues <= 0), 0)
  npc        <- length(evalues[evalues>0])
  npc        <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
  efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  # Estimated covariance function
  cov        <- efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  # Numerical integration for estimation of sigma2
  T.len      <- argvals[D] - argvals[1]  # total interval length
  T1.min     <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
  T1.max     <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
  DIAG       <- (diag.G0 - diag(cov))[T1.min:T1.max]  # function values
  w2         <- quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2     <- max(stats::weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
  ##
  
  #'# PACE ################################################################################
  if(PACE | PACE_E0){
    scoresP           <- vector(mode = "list", length(reconst_fcts))
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    wP                <- quadWeights(argvals, method = "trapezoidal")
    WsqrtP            <- diag(sqrt(wP))
    WinvsqrtP         <- diag(1/(sqrt(wP)))
    # Cov
    VP                <- WsqrtP %*% cov %*% WsqrtP
    evalP             <- eigen(VP, symmetric = TRUE, only.values = TRUE)$values
    evalP             <- replace(evalP, which(evalP <= 0), 0)
    npcP              <- length(evalP[evalP>0])  
    npcP              <- ifelse(is.null(pev), npcP, which(cumsum(evalP[evalP>0])/sum(evalP[evalP>0])>=pev)[1])
    efunctionsP       <- matrix(WinvsqrtP %*% eigen(VP, symmetric = TRUE)$vectors[, seq(len = npcP)], nrow = nrow(VP), ncol = npcP)
    evaluesP          <- evalP[1:npcP] 
    ##
    D.invP            <- diag(1/evaluesP, nrow = npcP, ncol = npcP)
    ##
    for(i in seq_len(length(reconst_fcts))){# i <- 1
      obs_locP        <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvals))
      Y.centP         <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
      ##
      if(sigma2 == 0){sigma2 <- 1e-6}
      if(PACE_E0){sigma2  <- 0}
      if(length(obs_locP) < npcP){npcP <- length(obs_locP)}
      ZcurP           <- efunctionsP[obs_locP, 1:npcP, drop=FALSE]
      ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # while(is.error(ZtZ_sD.invP)){# for preventing singularities
      #   sigma2          <- sigma2 + .Machine$double.eps
      #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # }
      if(is.error(ZtZ_sD.invP)){
        ZtZ_sD.invP     <- matrix(NA, npcP, npcP)
      }
      scoresP[[i]]    <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
    }
  } else {
    efunctionsP <- NA
    evaluesP    <- NA
    scoresP     <- NA  
  }
  #'# End PACE ############################################################################
  
  
  ## computations for observed fragments
  muO          <- vector("list", length(reconst_fcts))
  scoresO      <- vector("list", length(reconst_fcts))
  CE_scoresO    <- vector("list", length(reconst_fcts))
  evaluesO     <- vector("list", length(reconst_fcts))
  efunctionsO  <- vector("list", length(reconst_fcts))
  efun_reconst <- vector("list", length(reconst_fcts))
  ##
  obs_argvalsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){# i <- 1
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    w                 <- quadWeights(argvalsO[[i]], method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    locO              <- match(argvalsO[[i]],argvals)
    # CovOO
    VO                <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalO             <- eigen(VO, symmetric = TRUE, only.values = TRUE)$values
    evalO             <- replace(evalO, which(evalO <= 0), 0)
    npcO              <- length(evalO[evalO>0])  
    npcO              <- ifelse(is.null(pev), npcO, which(cumsum(evalO[evalO>0])/sum(evalO[evalO>0])>=pev)[1])
    efunctionsO[[i]]  <- matrix(Winvsqrt %*% eigen(VO, symmetric = TRUE)$vectors[, seq(len = npcO)], nrow = nrow(VO), ncol = npcO)
    evaluesO[[i]]     <- evalO[1:npcO]  
    ##
    D.inv             <- diag(1/evaluesO[[i]], nrow = npcO, ncol = npcO)
    Z                 <- efunctionsO[[i]]
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
    obs_locO          <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    obs_argvalsO[[i]] <- argvalsO[[i]][obs_locO]
    ## 
    if(CEscores){
      ## CEScores (i.e., PACE-Scores)
      if(sigma2           ==   0){sigma2 <- 1e-6}
      if(length(obs_locO) < npcO){npcO <- length(obs_locO)}
      ##
      Zcur           <- Z[obs_locO,1:npcO,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npcO,1:npcO])
      CE_scoresO[[i]] <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    } else {
      CE_scoresO[[i]] <- NA
    }
    ## Classical scores (via intergral approximation)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,,drop=FALSE], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=obs_argvalsO[[i]])})
    ##
    muO[[i]]     <- mu[locO]
    #'# ##################################################################
    ## Reconstructive eigenfunctions
    efun_reconst[[i]]  <- matrix(NA, nrow=length(argvals), ncol=npcO)
    ##
    for(k in seq_len(npcO)){
      efun_reconst[[i]][,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                     FUN = function(x){pracma::trapz(x=argvalsO[[i]], efunctionsO[[i]][,k] * x)})
      efun_reconst[[i]][,k] <- efun_reconst[[i]][,k] / evaluesO[[i]][k]
    }
    #'# ##################################################################
  }
  
  #'# Return results ##################################################
  ret.objects <- c("Y", "mu", "muO", "cov", "sigma2",
                   "argvals",    "argvalsO", "obs_argvalsO",
                   "CE_scoresO",  "scoresO", "scoresP",
                   "efunctions", "efunctionsO", "efun_reconst", "efunctionsP", 
                   "evalues",    "evaluesO", "evaluesP")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}

# AUXILIARY Functions -----------------------------------------------------

irreg2mat <- function(ydata, binning = FALSE, maxbins = 1000){
  ##
  ydata <- ydata[stats::complete.cases(ydata), ]
  nobs  <- length(unique(ydata$.id))
  newid <- as.numeric(as.factor(ydata$.id))
  bins  <- sort(unique(ydata$.index))
  if(binning && (length(bins) > maxbins)){
    binvalues <- seq((1 - 0.001 * sign(bins[1])) * bins[1], 
                     (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)], 
                     l = maxbins + 1)
    bins      <- binvalues
    binvalues <- utils::head(stats::filter(binvalues, c(0.5, 0.5)), -1)
  }else{
    binvalues <- bins
    bins      <- c((1 - 0.001 * sign(bins[1])) * bins[1], bins[-length(bins)], 
                   (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)])
    if(bins[1] == 0           ){bins[1]            <- -0.001}
    if(bins[length(bins)] == 0){bins[length(bins)] <-  0.001}
  }
  newindex         <- cut(ydata$.index, breaks = bins, include.lowest = TRUE)
  Y                <- matrix(NA, nrow = nobs, ncol = nlevels(newindex))
  colnames(Y)      <- binvalues
  attr(Y, "index") <- binvalues
  ## If there are more than one data-point within a bin, 
  ## then only one of these is used (the last one).
  Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
  ##
  return(Y)
}

quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method, 
                trapezoidal = {D <- length(argvals); 1/2 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)], argvals[D] - argvals[D - 1])}, 
                midpoint    = c(0, diff(argvals)), 
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule")
  )
  ##
  return(ret)
}

winsorize_x <- function(x, cut = 0.005){
  x <- x[!is.na(x)]
  ##
  cut_point_top    <- stats::quantile(x, 1 - cut)
  cut_point_bottom <- stats::quantile(x,     cut)
  i <-  which(x >= cut_point_top)
  j <-  which(x <= cut_point_bottom)
  x[i] <- cut_point_top
  x[j] <- cut_point_bottom
  ##
  return(x)
}


is.error <- function(x){
  bool.result <- inherits(x, "try-error")
  bool.result <- unname(bool.result)
  return(bool.result)
}


# Extrapolation ----------------------------------------------------------------
#' Reconstruct functional data through extrapolation
#'
#' This function allows to reconstruct partially observed functional data through an 
#' extrapolation. Each incomplete functional datum, observed up to a sample point T*,
#' is continuously extended with a straight line with slope the mean slope that all other
#' complete observations have from T* up to the right extrema of the complete observation
#' interval.
#' @param  curves   pxn matrix, p: number of discretization points, n=number of functions
#' @export out      a list containing the original and the reconstructed curves


extrapolation <- function(curves, t.points, T_hp, reconst_fcts)
{
  
  N <- length(t.points)
  ## Evaluation of mean slope
  mean.slope <- numeric(N-1)
  for(j in 1:(N-1))
  {
    mean.slope[j] <- mean(curves[N,-reconst_fcts]-curves[j,-reconst_fcts])/(t.points[N]-t.points[j])
  }
  
  ## All observations
  curves.rec <- curves
  
  for(i in 1:length(reconst_fcts))
  {
    Y         <- curves[,reconst_fcts[i]]
    id.na     <- is.na(Y)
    t.lastobs <- tail(t.points[!id.na], n=1)
    slope     <- tail(mean.slope[!id.na], n=1)
    Y.lastobs <- tail(Y[!id.na], n=1)
    curves.rec[id.na,reconst_fcts[i]] <- Y.lastobs + slope*(t.points[id.na]-t.lastobs)
  }
  
  ## Output
  
  out <- list(curves, curves.rec)
  names(out) <- c('curves', 'curves.rec')
  
  return(out)
  
}

