#' Conducting Cross-validated Kernel Ensemble based on user-specified formula.
#'
#' Perform Cross-validated Kernel Ensemble and optionally testing for kernel effect
#' based on user-specified formula.
#'
#' @param formula (formula) A user-supplied formula.
#' @param kern_func_list (list) a list of kernel functions in the kernel library
#' @param data (data.frame, n*d) a data.frame, list or environment (or object
#' coercible by as.data.frame to a data.frame), containing the variables in
#' formula. Neither a matrix nor an array will be accepted.
#' @param formula_test (formula) A user-supplied formula indicating the alternative
#' effect to test. All terms in the alternative mode must be specified as kernel terms.
#' @param mode (character) A character string indicating which tuning parameter
#' criteria is to be used.
#' @param strategy (character) A character string indicating which ensemble
#' strategy is to be used.
#' @param beta_exp (numeric/character) A numeric value specifying the parameter
#' when strategy = "exp" \code{\link{ensemble_exp}}.
#' @param lambda (numeric) A numeric string specifying the range of noise to be
#' chosen. The lower limit of lambda must be above 0.
#' @param verbose (logical) Whether to print additional messages.
#' @param test (character) Type of hypothesis test to conduct.
#' Must be eitehr 'asymp' or 'boot'.
#' @param alt_kernel_type (character) Type of alternative kernel effect to consider.
#' Must be either 'linear' or 'ensemble'
#' @param B (numeric) Number of boostrap samples.
#'
#' @return A list of results from kernel ensemble estimates. It contains
#' lambda (numeric) The estimated tuning parameter.
#' beta (matrix, d_fixed*1) fixed effect estimates.
#' alpha (matrix, n*1) kernel effect estimates.
#' K (matrix, n*n) Estimated ensemble kernel matrix.
#' u_hat (vector of length K) ensemble weight for kernel matrix.
#' base_est (list) result list output by the ensemble function.
#' @export cvek
#'
#' @examples
#' # create data
#' data <- as.data.frame(matrix(rnorm(700), ncol = 7,
#' dimnames = list(NULL, paste0("x", 1:7))))
#' data$y <- as.matrix(data) %*% rnorm(7)
#' formula <- y ~ x1 + x2 + k(x3, x4)
#' kern_par <- data.frame(method = c("rbf", "polynomial", "matern"),
#' l = c(.5, 1, 1.5), p = 1:3, stringsAsFactors = FALSE)
#'
#' # define kernel library
#' kern_func_list <- list()
#' for (d in 1:nrow(kern_par)) {
#' kern_func_list[[d]] <- generate_kernel(kern_par[d,]$method,
#' kern_par[d,]$l, kern_par[d,]$p)
#' }
#'
#' cvek(formula, kern_func_list = kern_func_list, data = data)
cvek <- function(formula,
                 kern_func_list,
                 data,
                 formula_test = NULL,
                 mode = "loocv",
                 strategy = "stack",
                 beta_exp = 1,
                 lambda = exp(seq(-10, 5)),
                 test = "boot",
                 alt_kernel_type = "linear",
                 B = 100,
                 verbose = FALSE) {
  # specify model matrices for main model
  model_matrices <- parse_cvek_formula(
    formula,
    kern_func_list = kern_func_list,
    data = data,
    verbose = verbose
  )
  
  # conduct estimation
  est_res <- estimation(
    Y = model_matrices$y,
    X = model_matrices$X,
    K_list = model_matrices$K,
    mode = mode,
    strategy = strategy,
    beta_exp = beta_exp,
    lambda = lambda
  )
  
  est_res$model_matrices <- model_matrices
  est_res$kern_func_list <- kern_func_list
  est_res$formula <- formula
  est_res$data <- data
  
  # conduct hypothesis test if formula_test is given.
  if (class(formula_test) == "formula") {
    est_res$pvalue <- cvek_test(est_res,
                                formula_test, 
                                kern_func_list,
                                data,
                                test = test,
                                alt_kernel_type = alt_kernel_type,
                                B = B,
                                verbose = verbose
    )
  }
  
  class(est_res) <- "cvek"
  est_res
}


#' Conduct Hypothesis Testing based on CVEK estimation result.
#'
#' @param est_res (list) Estimation results returned by estimation() procedure.
#' @param formula_test (formula) A user-supplied formula indicating the alternative
#' effect to test. All terms in the alternative mode must be specified as kernel terms.
#' @param kern_func_list (list) a list of kernel functions in the kernel library
#' @param data (data.frame, n*d) a data.frame, list or environment (or object
#' coercible by as.data.frame to a data.frame), containing the variables in
#' formula. Neither a matrix nor an array will be accepted.
#' @param test (character) Type of hypothesis test to conduct.
#' Must be eitehr 'asymp' or 'boot'.
#' @param alt_kernel_type (character) Type of alternative kernel effect to consider.
#' Must be either 'linear' or 'ensemble'
#' @param B (integer) A numeric value indicating times of resampling when test
#' = "boot".
#' @param verbose (logical) Whether to print additional messages.
#'
#' @return \item{pvalue}{(numeric) p-value of the test.}
#' @export
#'
#' @examples
#' @keywords internal
cvek_test <- function(est_res,
                      formula_test,
                      kern_func_list,
                      data,
                      test = "boot",
                      alt_kernel_type = "linear",
                      B = 100,
                      verbose = FALSE) {
  # TODO: generalize to allow ensemble effect kernels.
  alt_kernel_type <-
    match.arg(alt_kernel_type, c("linear", "ensemble"))
  
  # define kernel functions for alternative effect
  if (alt_kernel_type == "linear") {
    lnr_kern_func <- generate_kernel(method = "linear")
    alt_kern_func_list <- list(lnr_kern_func)
  } else if (alt_kernel_type == "ensemble") {
    alt_kern_func_list <- kern_func_list
  }
  
  # parse alternative effect formula to get kernel matrices
  test_matrices <- parse_cvek_formula(
    formula_test,
    kern_func_list = alt_kern_func_list,
    data = data,
    verbose = verbose
  )
  
  # check if alternative effect contain linear terms
  linear_terms <-
    setdiff(colnames(test_matrices$X), "(Intercept)")
  if (length(linear_terms) > 0) {
    stop(
      gettextf(
        "'formula_test' should contain only kernel terms. Linear terms found: %s",
        paste0(linear_terms, collapse = ", ")
      )
    )
  }
  
  # compute alternative kernel effect
  if (alt_kernel_type == "linear") {
    K_std_list <-
      lapply(test_matrices$K[[1]], function(K)
        K / tr(K))
    K_int <- Reduce("+", K_std_list)
  } else {
    u_weight <- est_res$u_hat
    K_int <- 0
    for (k in seq(length(kern_func_list))) {
      K_std_list <-
        lapply(test_matrices$K[[k]], function(K)
          K / tr(K))
      K_temp <- Reduce("+", K_std_list)
      K_int <- K_int + u_weight[k] * K_temp
    }
  }
  model_matrices <- est_res$model_matrices
  # estimate variance component parameters
  y_fixed <- model_matrices$X %*% est_res$beta
  sigma2_hat <- estimate_sigma2(
    Y = model_matrices$y,
    X = model_matrices$X,
    lambda_hat = est_res$lambda,
    y_fixed_hat = y_fixed,
    alpha_hat = est_res$alpha,
    K_hat = est_res$K
  )
  tau_hat <- sigma2_hat / est_res$lambda
  
  # compute p-value
  func_name <- paste0("test_", test)
  
  do.call(
    func_name, list(Y = model_matrices$y,
                    X = model_matrices$X,
                    K_int = K_int,
                    y_fixed = y_fixed,
                    alpha0 = est_res$alpha,
                    K_ens = est_res$K,
                    sigma2_hat = sigma2_hat,
                    tau_hat = tau_hat,
                    B = B
    )
  )
}

#' Parses user-supplied formula to fixed-effect and kernel matrices.
#'
#' @param formula (formula) A user-supplied formula.
#' @param kern_func_list (list) a list of kernel functions in the kernel library
#' @param data (data.frame, n*d) a data.frame, list or environment (or object
#' coercible by as.data.frame to a data.frame), containing the variables in
#' formula. Neither a matrix nor an array will be accepted.
#' @param data_new (data.frame, n_new*d) new data for computing predictions.
#' @param verbose (logical) Whether to print additional messages.
#'
#' @return A list of three slots:
#' \item{Y}{(matrix, n*1) The vector of response variable.}
#' \item{X}{(matrix, n*d_fix) The fixed effect matrix.}
#' \item{K}{(list of matrices) A nested list of
#' kernel term matrices. The first level corresponds to each base kernel
#' function in kern_func_list, the second level corresponds to each kernel term
#' specified in the formula.
#' }
#'
#' @details
#' The formula object is exactly like the formula for a GLM except that user can
#' use k() to specify kernel terms.
#' Additionally, user can specify interaction between kernel terms (using either '*' and ':'),
#' and exclude interaction term by including -1 on the rhs of formula.
#'
#' @author Jeremiah Zhe Liu
#' @export parse_cvek_formula
#'
#' @examples
#'
#' # create data
#' data <- as.data.frame(matrix(rnorm(700), ncol = 7,
#' dimnames = list(NULL, paste0("x", 1:7))))
#' data_new <- as.data.frame(matrix(rexp(700), ncol = 7,
#' dimnames = list(NULL, paste0("x", 1:7))))
#' data$y <- as.matrix(data) %*% rnorm(7)
#' formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6) + k(x7) + k(x3, x4) * k(x5, x6) * x7
#' kern_par <- data.frame(method = c("rbf", "polynomial", "matern"),
#' l = c(.5, 1, 1.5), d = 1:3, stringsAsFactors = FALSE)
#'
#' # define kernel library
#' kern_func_list <- list()
#' for (d in 1:nrow(kern_par)) {
#' kern_func_list[[d]] <- generate_kernel(kern_par[d,]$method,
#' kern_par[d,]$l, kern_par[d,]$d)
#' }
#'
#' # produce training data 
#' parse_cvek_formula(formula, kern_func_list = kern_func_list, 
#' data = data, data_new = NULL)
#' 
#' # produce prediction data from data_new
#' parse_cvek_formula(formula, kern_func_list = kern_func_list, 
#' data = data, data_new = data_new)
#'
parse_cvek_formula <-
  function(formula, kern_func_list, data, data_new = NULL, verbose = FALSE) {
    # extract dependent variables and terms
    tf <- terms.formula(formula, specials = c("k"))
    term_names <- attr(tf, "term.labels")
    num_terms <- length(term_names)
    
    var_table <- attr(tf, "factors")
    
    # extract indep variables and intercept
    intercept <- attr(tf, "intercept")
    response_vector <- NULL
    if ((attr(tf, "response") > 0) & is.null(data_new)) {
      response_name <- as.character(attr(tf, "variables")[2])
      response_vector <- data[, response_name]
    }
    
    # identify fixed-effect and kernel-effect terms
    kern_var_idx <- attr(tf, "specials")$k
    kern_term_idx <- NULL
    fixd_term_idx <- NULL
    
    if (length(kern_var_idx) > 0) {
      # if fomula contain kernel terms, identify their location
      kern_term_idx <-
        which(colSums(var_table[kern_var_idx, , drop = FALSE]) > 0)
      fixd_term_idx <- setdiff(1:num_terms, kern_term_idx)
      if (length(fixd_term_idx) == 0) {
        # set fixd_term_idx back to NULL if it is an empty set
        fixd_term_idx <- NULL
      }
    } else {
      fixd_term_idx <- 1:num_terms
    }
    
    # assemble fixed-effect and kernel-effect formula
    kern_effect_formula <- NULL
    fixed_effect_formula <- NULL
    
    if (!is.null(kern_term_idx)) {
      kern_effect_formula <-
        as.formula(paste("~", paste(term_names[kern_term_idx], collapse = " + ")))
    }
    
    if (!is.null(fixd_term_idx)) {
      fixed_effect_formula <-
        paste("~", paste(term_names[fixd_term_idx], collapse = " + "))
      if (intercept == 0)
        fixed_effect_formula <- paste0(fixed_effect_formula, " -1")
      
      fixed_effect_formula <- as.formula(fixed_effect_formula)
    } else if (intercept > 0) {
      # intercept only
      fixed_effect_formula <- ~ 1
    }
    
    # produce fixed-effect matrix X using model.frame
    fixed_effect_matrix <- NULL
    if (!is.null(fixed_effect_formula)) {
      fixed_effect_data <- data
      if (!is.null(data_new)){
        fixed_effect_data <- data_new
      }
      
      fixed_effect_matrix <-
        model.matrix(fixed_effect_formula, data = fixed_effect_data)
    }
    
    # produce kernel-effect matrices Ks
    # (list of kernel terms, one for each kern in library)
    kernel_effect_matrix_list <-
      vector("list", length = length(kern_func_list))
    names(kernel_effect_matrix_list) <- names(kern_func_list)
    
    if (!is.null(kern_effect_formula)) {
      if (verbose)
        print("Preparing Kernels...")
      for (kern_func_id in 1:length(kern_func_list)) {
        # prepare kernel function
        kern_func <- kern_func_list[[kern_func_id]]
        
        # compute kernel terms using the kern_func, then store to list
        kernel_effect_matrix_list[[kern_func_id]] <-
          parse_kernel_terms(kern_effect_formula, kern_func, 
                             data = data, data_new = data_new)
      }
      if (verbose)
        print("Done!")
    }
    
    # combine and return
    # TODO (jereliu): discuss whether to add standardization
    list(y = response_vector,
         X = fixed_effect_matrix,
         K = kernel_effect_matrix_list)
  }

#' Computes kernel matrix for each kernel term in the formula
#'
#' @param kern_effect_formula (character) a term in the formula
#' @param kern_func (function) A kernel function. Will be overwritten to linear
#' kernel if the variable doesn't contain 'k()'
#' @param data (data.frame, n*d) a data.frame, list or environment (or object
#' coercible by as.data.frame to a data.frame), containing the variables in
#' formula. Neither a matrix nor an array will be accepted.
#' @param data_new (data.frame, n_new*d) new data for computing predictions.

#' @return A list of kernel matrices for each term in the formula
#' @export parse_kernel_terms
#'
#' @author Jeremiah Zhe Liu
#' @keywords internal
parse_kernel_terms <-
  function(kern_effect_formula, kern_func, data, data_new = NULL) {
    kern_var_table <- attr(terms(kern_effect_formula), "factors")
    kern_var_names <- rownames(kern_var_table)
    kern_term_names <- colnames(kern_var_table)
    
    # prepare kernel variables
    kern_var_list <- vector("list", length = length(kern_var_names))
    names(kern_var_list) <- kern_var_names
    
    for (var_name in kern_var_names) {
      # compute individual kernel variables for each base kernel
      kern_var_list[[var_name]] <-
        parse_kernel_variable(var_name, kern_func = kern_func, 
                              data = data, data_new = data_new)
    }
    
    # assemble kernel terms (for interaction terms)
    kern_term_list <-
      vector("list", length = length(kern_term_names))
    names(kern_term_list) <- kern_term_names
    
    for (term_id in 1:length(kern_term_names)) {
      term_name <- kern_term_names[term_id]
      kern_term_var_idx <- which(kern_var_table[, term_id] > 0)
      # compute interaction term
      kern_term_list[[term_name]] <-
        Reduce("*", kern_var_list[kern_term_var_idx])
    }
    
    kern_term_list
  }

#' Creates kernel matrix for each variable in the formula
#'
#' @param kern_var_name (vector of characters) Names of variables in data to
#' create the kernel matrix from. Must be a single term that is either of the form
#' \eqn{x} (a single linear term) or \eqn{k(x1, x2, \dots)} (a kernel term that may
#' contain multiple variables).
#' @param kern_func (function) A kernel function. Will be overwritten to linear
#' kernel if the variable doesn't contain 'k()'
#' @param data (data.frame, n*d) a data.frame, list or environment (or object
#' coercible by as.data.frame to a data.frame), containing the variables in
#' formula. Neither a matrix nor an array will be accepted.
#' @param data_new (data.frame, n_new*d) new data for computing predictions.
#'
#' @return The n*n kernel matrix corresponding to the variable being computed.
#' @export parse_kernel_variable
#'
#' @author Jeremiah Zhe Liu
#' @keywords internal
parse_kernel_variable <- function(kern_var_name, kern_func, data, data_new=NULL) {
  # parse kernel term
  kern_term_formula <-
    terms.formula(formula(paste("~", kern_var_name)),
                  specials = c("k"))
  
  is_fixed_eff <- is.null(attr(kern_term_formula, "specials")$k)
  input_var_names <- all.vars(kern_term_formula)
  
  # sets kern_func to linear kernel if kern_var is a fixed-effect term
  if (is_fixed_eff) {
    kern_func <- generate_kernel(method = "linear")
  }
  
  # extract data
  Z_extract_command <- gettextf("cbind(%s)",
                                paste(input_var_names, collapse = ", "))
  Z_mat <- eval(parse(text = Z_extract_command), envir = data)
  
  Z_mat_new <- Z_mat
  
  if (!is.null(data_new)){
    Z_mat_new <- eval(parse(text = Z_extract_command), envir = data_new)
  }
  
  # compute kernel matrix and return
  kern_func(Z_mat_new, Z_mat)
}