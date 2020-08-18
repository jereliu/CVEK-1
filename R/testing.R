#' Conducting Score Tests for Interaction Using Asymptotic Test
#' 
#' Conduct score tests comparing a fitted model and a more general alternative
#' model using asymptotic test.
#' 
#' \bold{Asymptotic Test}
#' 
#' This is based on the classical variance component test to construct a
#' testing procedure for the hypothesis about Gaussian process function.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param y_fixed (vector of length n) Estimated fixed effect of the
#' response.
#' @param alpha0 (vector of length n) Kernel effect estimator of the estimated
#' ensemble kernel matrix.
#' @param K_ens (matrix, n*n) Estimated ensemble kernel matrix.
#' @param K_int (matrix, n*n) The kernel matrix to be tested.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effect.
#' @param tau_hat (numeric) The estimated noise of the kernel effect.
#' @param B (integer) A numeric value indicating times of resampling when test
#' = "boot".
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#' 
#' mode: \code{\link{tuning}}
#' 
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#' 
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#' 
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
test_asymp <- function(Y, X, y_fixed, alpha0, K_ens, 
                      K_int, sigma2_hat, tau_hat, B) {
  
  n <- length(Y)
  score_chi <-
    compute_stat(Y, K_int, y_fixed, K_ens, sigma2_hat, tau_hat)
  K0 <- K_ens
  V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
  P0_mat <- V0_inv - V0_inv %*%
    X %*% ginv(t(X) %*% V0_inv %*% X) %*% t(X) %*% V0_inv
  drV0_tau <- K0
  drV0_sigma2 <- diag(n)
  drV0_del <- tau_hat * K_int
  I0 <- compute_info(P0_mat,
                     mat_del = drV0_del, mat_sigma2 = drV0_sigma2,
                     mat_tau = drV0_tau)
  tot_dim <- ncol(I0)
  I_deldel <-
    I0[1, 1] -
    I0[1, 2:tot_dim] %*% ginv(I0[2:tot_dim, 2:tot_dim]) %*% I0[2:tot_dim, 1]
  md <- tau_hat * sum(diag(K_int %*% P0_mat)) / 2
  m_chi <- I_deldel / (2 * md)
  d_chi <- md / m_chi
  pvalue <- 1 - pchisq(score_chi / m_chi, d_chi)
  
  pvalue
}





#' Conducting Score Tests for Interaction Using Bootstrap Test
#' 
#' Conduct score tests comparing a fitted model and a more general alternative
#' model using bootstrap test.
#' 
#' \bold{Bootstrap Test}
#' 
#' When it comes to small sample size, we can use bootstrap test instead, which
#' can give valid tests with moderate sample sizes and requires similar
#' computational effort to a permutation test.
#' 
#' @param Y (matrix, n*1) The vector of response variable.
#' @param X (matrix, n*d_fix) The fixed effect matrix.
#' @param y_fixed (vector of length n) Estimated fixed effect of the
#' response.
#' @param alpha0 (vector of length n) Kernel effect estimator of the estimated
#' ensemble kernel matrix.
#' @param K_ens (matrix, n*n) Estimated ensemble kernel matrix.
#' @param K_int (matrix, n*n) The kernel matrix to be tested.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effect.
#' @param tau_hat (numeric) The estimated noise of the kernel effect.
#' @param B (integer) A numeric value indicating times of resampling when test
#' = "boot".
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{generate_kernel}}
#' 
#' mode: \code{\link{tuning}}
#' 
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#' 
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#' 
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
test_boot <- function(Y, X, y_fixed, alpha0, K_ens, 
                      K_int, sigma2_hat, tau_hat, B) {
  
  n <- length(Y)
  meanY <- K_ens %*% alpha0 + y_fixed
  bs_test <- sapply(1:B, function(k) {
    Ystar <- meanY + rnorm(n, sd = sqrt(sigma2_hat))
    compute_stat(Ystar, K_int, y_fixed, K_ens, sigma2_hat, tau_hat)
  })
  original_test <-
    compute_stat(Y, K_int, y_fixed, K_ens, sigma2_hat, tau_hat)
  pvalue <- mean(as.numeric(original_test) <= bs_test)
  
  pvalue
}
