#' Conducting Gaussian Process Regression
#' 
#' Conduct Gaussian process regression based on the estimated ensemble kernel
#' matrix.
#' 
#' After obtaining the ensemble kernel matrix, we can calculate the output of
#' Gaussian process regression.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param K_list (list of matrices) A nested list of kernel term matrices. 
#' The first level corresponds to each base kernel function in kern_func_list, 
#' the second level corresponds to each kernel term specified in the formula.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param strategy (character) A character string indicating which ensemble
#' strategy is to be used.
#' @param beta_exp (numeric/character) A numeric value specifying the parameter
#' when strategy = "exp" \code{\link{ensemble_exp}}.
#' @param lambda (numeric) A numeric string specifying the range of tuning 
#' parameter to be chosen. The lower limit of lambda must be above 0.
#' @param ... Additional parameters to pass to estimate_ridge.
#' 
#' @return \item{lambda}{(numeric) The selected tuning parameter based on the
#' estimated ensemble kernel matrix.}
#' 
#' \item{beta}{(matrix, d_fixed*1) Fixed effect estimates.}
#' 
#' \item{alpha}{(matrix, n*1) Kernel effect estimates.}
#' 
#' \item{K}{(matrix, n*n) Estimated ensemble kernel matrix.}
#' 
#' \item{u_hat}{(vector of length K) A vector of weights of the kernels in the
#' library.}
#' 
#' \item{kern_term_effect}{(matrix, n*n) Estimated ensemble kernel effect matrix.}
#' 
#' \item{base_est}{(list) The detailed estimation results of K kernels.}
#'  
#' @author Wenying Deng
#' @seealso strategy: \code{\link{ensemble}}
#' 
#' @export estimation
estimation <- function(Y, X, K_list = NULL,
                       mode = "loocv", 
                       strategy = "stack", 
                       beta_exp = 1,
                       lambda = exp(seq(-10, 5)), 
                       ...) {
  
  # base model estimate
  n <- length(Y)
  if(sum(X[, 1] == 1) != n) {
    X <- cbind(matrix(1, nrow = n, ncol = 1), X)
  }
  
  kern_size <- length(K_list)
  if (kern_size == 0) {
    # if K_list is empty, estimate only fixed effect
    # 'assemble' fake ensemble kernel matrix
    base_est <- NULL
    K_ens <- NULL
    
    # final estimate
    lambda_ens <- tuning(Y, X, K_ens, mode, lambda)
    ens_est <- estimate_ridge(Y = Y, X = X, K = K_ens, lambda = lambda_ens, ...)
    
    # kernel terms estimates
    u_weight <- 1
    kern_term_effect <- NULL
  } else {
    # if K_list is not empty, estimate full ensemble 
    base_est <- estimate_base(Y, X, K_list, mode, lambda, ...)
    P_K_hat <- base_est$P_K_hat
    error_mat <- base_est$error_mat
    ens_res <- ensemble(strategy, beta_exp, error_mat, P_K_hat)
    
    # assemble ensemble kernel matrix
    K_ens <- ensemble_kernel_matrix(ens_res$A_est)
    K_ens <- list(K_ens)
    
    # final estimate
    lambda_ens <- tuning(Y, X, K_ens, mode, lambda)
    ens_est <- estimate_ridge(Y = Y, X = X, K = K_ens, lambda = lambda_ens, ...)
    
    # kernel terms estimates
    u_weight <- ens_res$u_hat
    kern_term_effect <- 0
    for (k in seq(kern_size)) {
      kern_term_effect <- kern_term_effect + 
        u_weight[k] * base_est$kern_term_list[[k]]
    }
  }
  
  list(lambda = lambda_ens, beta = ens_est$beta, 
       alpha = ens_est$alpha, K = K_ens[[1]], 
       u_hat = u_weight, 
       kern_term_effect = kern_term_effect, 
       base_est = base_est)
}



#' Estimating Projection Matrices
#' 
#' Calculate the estimated projection matrices for every kernels in the kernel
#' library.
#' 
#' For a given mode, this function returns a list of projection matrices for
#' every kernel in the kernel library and a n*K matrix indicating errors.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param K_list (list of matrices) A nested list of kernel term matrices. 
#' The first level corresponds to each base kernel function in kern_func_list, 
#' the second level corresponds to each kernel term specified in the formula.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param lambda (numeric) A numeric string specifying the range of tuning 
#' parameter to be chosen. The lower limit of lambda must be above 0.
#' @param ... Additional parameters to pass to estimate_ridge.
#' 
#' @return \item{A_hat}{(list of length K) A list of projection matrices for
#' each kernel in the kernel library.}
#' 
#' \item{P_K_hat}{(list of length K) A list of projection matrices to kernel
#' space for each kernel in the kernel library.}
#' 
#' \item{beta_list}{(list of length K) A list of fixed effect estimators for
#' each kernel in the kernel library.}
#' 
#' \item{alpha_list}{(list of length K) A list of kernel effect estimates for
#' each kernel in the kernel library.}
#' 
#' \item{kern_term_list}{(list of length K) A list of kernel effects for
#' each kernel in the kernel library.}
#' 
#' \item{A_proc_list}{(list of length K) A list of projection matrices for
#' each kernel in the kernel library.}
#' 
#' \item{lambda_list}{(list of length K) A list of selected tuning parameters for
#' each kernel in the kernel library.}
#' 
#' \item{error_mat}{(matrix, n*K) A n\*K matrix indicating errors.}
#' 
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaussian Processes. October 2017.
#' @keywords internal
#' @export estimate_base
estimate_base <- function(Y, X, K_list, mode, lambda, ...) {
  A_hat <- list()
  P_K_hat <- list()
  beta_list <- list()
  alpha_list <- list()
  lambda_list <- list()
  kern_term_list <- list()
  A_proc_list <- list()
  n <- length(Y)
  kern_size <- length(K_list)
  error_mat <- matrix(0, nrow = n, ncol = kern_size)
  
  for (k in seq(kern_size)) {
    lambda0 <- tuning(Y, X, K_list[[k]], mode, lambda)
    estimate <- estimate_ridge(Y = Y, X = X, 
                               K = K_list[[k]], 
                               lambda = lambda0, ...)
    A <- estimate$proj_matrix$total
    
    # produce loocv error matrix
    error_mat[, k] <- (diag(n) - A) %*% Y / (1 - diag(A))
    
    A_hat[[k]] <- A  
    P_K_hat[[k]] <- estimate$proj_matrix$P_K0
    beta_list[[k]] <- estimate$beta
    alpha_list[[k]] <- estimate$alpha
    kern_term_list[[k]] <- estimate$kern_term_mat
    lambda_list[[k]] <- lambda0
    A_proc_list[[k]] <- estimate$A_list
  }
  
  list(A_hat = A_hat, P_K_hat = P_K_hat,
       beta_list = beta_list, alpha_list = alpha_list,
       kern_term_list = kern_term_list, 
       A_proc_list = A_proc_list, 
       lambda_list = lambda_list, error_mat = error_mat)
}



#' Estimating a Single Model
#' 
#' Estimating projection matrices and parameter estimates for a single model.
#' 
#' For a single model, we can calculate the output of gaussian process
#' regression, the solution is given by \deqn{\hat{\beta}=[X^T(K+\lambda
#' I)^{-1}X]^{-1}X^T(K+\lambda I)^{-1}y} \deqn{\hat{\alpha}=(K+\lambda
#' I)^{-1}(y-\hat{\beta}X)}.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param K (list of matrices) A nested list of kernel term matrices, 
#' corresponding to each kernel term specified in the formula for 
#' a base kernel function in kern_func_list.
#' @param lambda (numeric) A numeric string specifying the range of tuning parameter 
#' to be chosen. The lower limit of lambda must be above 0.
#' @param compute_kernel_terms (logic) Whether to computing effect for each individual terms.
#' If FALSE then only compute the overall effect.
#' @param converge_thres (numeric) The convergence threshold for computing kernel terms.
#' 
#' @return \item{beta}{(matrix, d_fixed*1) Fixed effect estimates.}
#' 
#' \item{alpha}{(matrix, n*k_terms) Kernel effect estimates for each kernel term.}
#' 
#' \item{kern_term_mat}{(matrix, n*k_terms) Kernel effect for each kernel term.}
#' 
#' \item{A_list}{(list of length k_terms) Projection matrices for each kernel term.}
#' 
#' \item{proj_matrix}{(list of length 4) Estimated projection matrices, combined 
#' across kernel terms.}
#' @author Wenying Deng
#' @references Andreas Buja, Trevor Hastie, and Robert Tibshirani. (1989) 
#' Linear Smoothers and Additive Models. Ann. Statist. Volume 17, Number 2, 453-510.
#' @examples
#' kern_par <- data.frame(method = c("linear", "polynomial", "matern"), 
#'                        l = c(.5, 1, 1.5), p = 1:3, stringsAsFactors = FALSE)
#' # define kernel library
#' kern_func_list <- define_library(kern_par)
#' n <- 50
#' d <- 4
#' formula <- y ~ x1 + x2 + k(x3, x4)
#' data <- as.data.frame(matrix(
#'   rnorm(n * d),
#'   ncol = d,
#'   dimnames = list(NULL, paste0("x", 1:d))
#' ))
#' lnr_kern_func <- generate_kernel(method = "linear")
#' kern_effect_mat <- 
#'   parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
#' beta_true <- c(1, .41, 2.37)
#' alpha_true <- rnorm(n)
#' 
#' K <- kern_effect_mat
#' data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
#'   K %*% alpha_true
#' 
#' model_matrices <- parse_cvek_formula(formula, 
#'                                      kern_func_list = kern_func_list, 
#'                                      data = data, verbose = FALSE)
#' estimate_ridge(Y = model_matrices$y, X = model_matrices$X, 
#' K = model_matrices$K[[1]], lambda = exp(-5))
#' 
#' @export estimate_ridge
estimate_ridge <- function(Y, X, K, lambda, 
                           compute_kernel_terms=TRUE, 
                           converge_thres=1e-4){
  # standardize kernel matrix
  n <- length(Y)
  if(sum(X[, 1] == 1) != n) {
    X <- cbind(matrix(1, nrow = n, ncol = 1), X)
  }
  
  # initialize parameters and calculate projection matrices
  X_mat <- ginv(t(X) %*% X) %*% t(X)
  beta <- X_mat %*% Y
  P_K <- 0
  P_X <- X %*% X_mat
  A_list <- NULL
  P_K0 <- NULL
  alpha_mat <- NULL
  kern_term_mat <- NULL
  
  if (!is.null(K)) {
    # prepare matrices needed for computing project
    alpha_mat <- matrix(0, nrow = n, ncol = length(K))
    A <- 0
    H <- X %*% X_mat
    V_inv_list <- list()
    A_list <- list()
    for (d in seq(length(K))) {
      K[[d]] <- K[[d]] / sum(diag(K[[d]]))
      V_inv_list[[d]] <- ginv(K[[d]] + lambda * diag(n))
      alpha_mat[, d] <- 
        V_inv_list[[d]] %*% (Y - X %*% beta)
      S_d <- K[[d]] %*% V_inv_list[[d]]
      A_list[[d]] <- (diag(n) + K[[d]] / lambda) %*% S_d
      A <- A + A_list[[d]]
    }
    B <- ginv(diag(n) + A) %*% A
    
    # project matrices
    # residual projection to kernel space
    P_K <- ginv(diag(n) - B %*% H) %*% B %*% (diag(n) - H)
    # projection to fixed-effect space
    P_X <- H %*% (diag(n) - P_K)
    # projection to kernel space
    # P_K0 <- P_K %*% ginv(diag(n) - P_X)
    P_K0 <- B
    
    # compute kernel term effects 
    if (!compute_kernel_terms) {
      # compute overall effect
      beta <- X_mat %*% (diag(n) - P_K) %*% Y
      kern_term_mat <- P_K %*% Y
    } else {
      # compute individual terms by iterative update using backfitting algorithm
      alpha_temp <- 1e3 * alpha_mat
      
      beta <- X_mat %*% (diag(n) - P_K) %*% Y
      while (euc_dist(alpha_mat, alpha_temp) > converge_thres) {
        alpha_temp <- alpha_mat
        kernel_effect <- 0
        for (d in seq(length(K))) {
          kernel_effect <- kernel_effect + K[[d]] %*% alpha_mat[, d]
        }
        for (d in seq(length(K))) {
          alpha_mat[, d] <- 
            V_inv_list[[d]] %*% 
            (Y - X %*% beta - kernel_effect + K[[d]] %*% alpha_mat[, d])
        }
      }      
      
      # kernel terms estimates
      kern_term_mat <- matrix(0, nrow = n, ncol = length(K))
      for (d in seq(length(K))) {
        kern_term_mat[, d] <- K[[d]] %*% alpha_mat[, d]
      }
    }
    
  }
  proj_matrix_list <- list(total = P_X + P_K,
                           P_X = P_X, P_K = P_K,
                           P_K0 = P_K0)
  
  list(beta = beta, alpha = alpha_mat, 
       kern_term_mat = kern_term_mat, 
       A_list = A_list, 
       proj_matrix = proj_matrix_list)
}




#' Calculating Ensemble Kernel Matrix
#' 
#' Calculating ensemble kernel matrix and truncating those columns whose
#' eigenvalues are smaller than the given threshold.
#' 
#' After we obtain the ensemble projection matrix, we can calculate the
#' ensemble kernel matrix.
#' 
#' @param A_est (matrix) Ensemble projection matrix.
#' @param eig_thres (numeric) Threshold to truncate the kernel matrix.
#' @return \item{K_hat}{(matrix, n*n) Estimated ensemble kernel matrix.}
#' @author Wenying Deng
#' @keywords internal
#' @export ensemble_kernel_matrix
ensemble_kernel_matrix <- function(A_est, eig_thres = 1e-11) {
  As <- svd(A_est)
  U <- As$u
  d <- As$d
  
  # produce spectral components for ensemble kernel matrix
  ensemble_dim <- sum(d > eig_thres)
  U_ens <- U[, 1:ensemble_dim]
  d_ens <- d[1:ensemble_dim] / (1 - d[1:ensemble_dim])
  
  # assemble ensemble matrix and return
  K_hat <- U_ens %*% diag(d_ens) %*% t(U_ens)
  K_hat / sum(diag(K_hat))
}
