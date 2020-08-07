#' Defining Kernel Library
#' 
#' Generate the expected kernel library based on user-specified dataframe.
#' 
#' It creates a kernel library according to the parameters given in kern_par.
#' 
#' * kern_par: for a library of K kernels, the dimension of this dataframe is
#'   K*3. Each row represents a kernel. The first column is method, with entries
#'   of character class. The second and the third are l and p respectively, both
#'   with entries of numeric class.
#' 
#' @param kern_par (dataframe, K*3) A dataframe indicating the parameters of
#' base kernels to fit kernel effect. See Details.
#' @return \item{kern_func_list}{(list of length K) A list of kernel functions 
#' given by user. Will be overwritten to linear kernel if kern_par is NULL.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#' 
#' @export define_library
define_library <- function(kern_par = NULL) {
  if (is.null(kern_par)){
    lnr_func <- generate_kernel(method = "linear")
    kern_func_list <- list(lnr_func)
  } else {
    kern_func_list <- list()
    for (d in 1:nrow(kern_par)) {
      kern_func_list[[d]] <- generate_kernel(kern_par[d,]$method,
                                             kern_par[d,]$l, 
                                             kern_par[d,]$p,
                                             kern_par[d,]$sigma)
    }
  }
  
  kern_func_list
}



#' Estimating Noise
#' 
#' An implementation of Gaussian processes for estimating noise.
#' 
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param lambda_hat (numeric) The selected tuning parameter based on the
#' estimated ensemble kernel matrix.
#' @param y_fixed_hat (vector of length n) Estimated fixed effect of the
#' response.
#' @param alpha_hat (vector of length n) Kernel effect estimators of the
#' estimated ensemble kernel matrix.
#' @param K_hat (matrix, n*n) Estimated ensemble kernel matrix.
#' @return \item{sigma2_hat}{(numeric) The estimated noise of the fixed
#' effect.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaussian Processes. October 2017.s
#' @keywords internal
#' @export estimate_sigma2
estimate_sigma2 <- function(Y, X, lambda_hat, y_fixed_hat, alpha_hat, K_hat) {

  n <- length(Y)
  V_inv <- ginv(K_hat + lambda_hat * diag(n))
  B_mat <- ginv(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
  P_X <- X %*% B_mat
  P_K <- K_hat %*% V_inv %*% (diag(n) - P_X)
  A <- P_X + P_K
  sigma2_hat <- sum((Y - y_fixed_hat - K_hat %*% alpha_hat) ^ 2) / (n - sum(diag(A)) - 1)
  
  sigma2_hat
}



#' Computing Score Test Statistics.
#' 
#' Compute score test statistics.
#' 
#' The test statistic is distributed as a scaled Chi-squared distribution.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param K_int (matrix, n*n) The kernel matrix to be tested.
#' @param y_fixed (vector of length n) Estimated fixed effect of the
#' response.
#' @param K_0 (matrix, n*n) Estimated ensemble kernel matrix.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effect.
#' @param tau_hat (numeric) The estimated noise of the kernel effect.
#' @return \item{test_stat}{(numeric) The computed test statistic.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @keywords internal
#' @export compute_stat
compute_stat <-
  function(Y, K_int, y_fixed, K0, sigma2_hat, tau_hat) {
    
    n <- length(Y)
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    test_stat <- tau_hat * t(Y - y_fixed) %*% V0_inv %*%
      K_int %*% V0_inv %*% (Y - y_fixed) / 2

    test_stat
  }



#' Computing Information Matrices
#' 
#' Compute information matrices based on block matrices.
#' 
#' This function gives the information value of the interaction strength.
#' 
#' @param P0_mat (matrix, n*n) Scale projection matrix under REML.
#' @param mat_del (matrix, n*n) Derivative of the scale covariance matrix of Y
#' with respect to delta.
#' @param mat_sigma2 (matrix, n*n) Derivative of the scale covariance matrix of
#' Y with respect to sigma2.
#' @param mat_tau (matrix, n*n) Derivative of the scale covariance matrix of Y
#' with respect to tau.
#' @return \item{I0}{(matrix, n*n) The computed information value.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @keywords internal
#' @export compute_info
compute_info <-
  function(P0_mat, mat_del = NULL, mat_sigma2 = NULL, mat_tau = NULL) {
    
    I0 <- matrix(NA, 3, 3)
    I0[1, 1] <- sum(diag(P0_mat %*% mat_del %*% P0_mat %*% mat_del)) / 2  
    I0[1, 2] <- sum(diag(P0_mat %*% mat_del %*% P0_mat %*% mat_sigma2)) / 2
    I0[2, 1] <- I0[1, 2]
    I0[1, 3] <- sum(diag(P0_mat %*% mat_del %*% P0_mat %*% mat_tau)) / 2
    I0[3, 1] <- I0[1, 3]
    I0[2, 2] <- sum(diag(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_sigma2)) / 2
    I0[2, 3] <- sum(diag(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_tau)) / 2  
    I0[3, 2] <- I0[2, 3]
    I0[3, 3] <- sum(diag(P0_mat %*% mat_tau %*% P0_mat %*% mat_tau)) / 2  

    I0
  }



#' Standardizing Matrix
#' 
#' Center and scale the data matrix into mean zero and standard deviation one.
#' 
#' This function gives the standardized data matrix.
#' 
#' @param X (matrix) Original data matrix.
#' @return \item{X}{(matrix) Standardized data matrix.}
#' @author Wenying Deng
#' @keywords internal
#' @export standardize
standardize <- function(X) {
  
  Xm <- colMeans(X)
  n <- nrow(X)
  p <- ncol(X)
  X <- X - rep(Xm, rep(n, p))
  Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
  X <- X / rep(Xscale, rep(n, p))
  
  X
}



#' Computing Euclidean Distance between Two Vectors (Matrices)
#' 
#' Compute the L2 distance between two vectors or matrices.
#' 
#' This function gives the Euclidean distance between two 
#' vectors or matrices.
#' 
#' @param x1 (vector/matrix) The first vector/matrix.
#' @param x2 (vector/matrix, the same dimension as x1) 
#' The second vector/matrix.
#' @return \item{dist}{(numeric) Euclidean distance.}
#' @author Wenying Deng
#' @keywords internal
#' @export euc_dist
euc_dist <- function(x1, x2 = NULL) {
  if (is.null(x2)) {
    dist <- sqrt(sum(x1 ^ 2))
  } else {
    dist <- sqrt(sum((x1 - x2) ^ 2))
  }
  dist
}
