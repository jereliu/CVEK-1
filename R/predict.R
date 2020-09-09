#' Predicting New Response
#' 
#' Predicting new response based on given design matrix and 
#' the estimation result.
#' 
#' After we obtain the estimation result, we can predict new response.
#' 
#' @param object (list) Estimation results returned by cvek() procedure.
#' @param newdata (dataframe) The new set of predictors, whose name is 
#' the same as those of formula in cvek().
#' @param ... Further arguments passed to or from other methods.
#' @return \item{y_pred}{(matrix, n*1) Predicted new response.}
#' @author Wenying Deng
#' @examples
#'
#' kern_par <- data.frame(method = rep("rbf", 3),
#' l = rep(3, 3), p = rep(2, 3), 
#' stringsAsFactors = FALSE)
#' # define kernel library
#' kern_func_list <- define_library(kern_par)
#' 
#' n <- 100
#' d <- 4
#' formula <- y ~ x1 + x2 + k(x3, x4)
#' formula_test <- y ~ k(x1, x2) * k(x3, x4)
#' data <- as.data.frame(matrix(
#'   rnorm(n * d),
#'   ncol = d,
#'   dimnames = list(NULL, paste0("x", 1:d))
#' ))
#' beta_true <- c(1, .41, 2.37)
#' lnr_kern_func <- generate_kernel(method = "rbf", l = 3)
#' kern_effect_lnr <-
#'   parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
#' alpha_lnr_true <- rnorm(n)
#' 
#' data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true +
#'   kern_effect_lnr %*% alpha_lnr_true
#' 
#' data_train <- data[1:50, ]
#' data_test <- data[51:100, ]
#' 
#' result <- cvek(formula,
#'                kern_func_list,
#'                data_train,
#'                formula_test,
#'                mode = "loocv",
#'                strategy = "stack",
#'                beta_exp = 1,
#'                lambda = exp(seq(-10, 5)),
#'                test = "boot",
#'                alt_kernel_type = "linear",
#'                B = 100,
#'                verbose = FALSE)
#' 
#' predict(result, data_test)
#' 
#' @importFrom utils data
#' @export predict.cvek
#' @export
predict.cvek <- function(object, newdata, ...) {
  
  model_matrices <- object$model_matrices
  kern_func_list <- object$kern_func_list
  new_matrices <- parse_cvek_formula(object$formula, 
                                     kern_func_list, 
                                     data = object$data, 
                                     data_new = newdata)
  X <- new_matrices$X
  A <- 0
  Xmat <- ginv(t(model_matrices$X) 
               %*% model_matrices$X) %*% t(model_matrices$X)
  H <- model_matrices$X %*% Xmat
  H_star <- X %*% Xmat
  n <- length(object$alpha)
  P_K_star <- list()
  P_X_star <- list()
  y_pred <- 0
  for (k in seq(length(kern_func_list))) {
    
    B_temp <- 0
    for (d in seq(length(new_matrices$K[[k]]))) {
      S_d_star <- new_matrices$K[[k]][[d]] %*% ginv(model_matrices$K[[k]][[d]] 
                                                    + object$base_est$lambda_list[[k]] * diag(n))
      B_temp <- B_temp + S_d_star %*% (diag(n) + object$base_est$A_proc_list[[k]][[d]])
    }
    B_star <- B_temp %*% (diag(n) - object$base_est$P_K_hat[[k]])
    P_K <- ginv(diag(n) - object$base_est$P_K_hat[[k]] %*% H) %*% 
      object$base_est$P_K_hat[[k]] %*% (diag(n) - H)
    P_K_star[[k]] <- B_star %*% (diag(n) - H + H %*% object$base_est$P_K_hat[[k]])
    P_X_star[[k]] <- H_star %*% (diag(n) - P_K)
    y_pred <- y_pred + object$u_hat[k] * (P_K_star[[k]] + P_X_star[[k]]) %*% object$model_matrices$y
  }

  y_pred
}
