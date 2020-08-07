context("Hypothesis testing")

test_that("range of output", {
  kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                         l = c(.5, 1, 1.5), p = 1:3, stringsAsFactors = FALSE)
  # define kernel library
  kern_func_list <- define_library(kern_par)
  n <- 50
  d <- 6
  formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6)
  data <- as.data.frame(matrix(
    rnorm(n * d),
    ncol = d,
    dimnames = list(NULL, paste0("x", 1:d))
  ))
  lnr_kern_func <- generate_kernel(method = "linear")
  rbf_kern_func <- generate_kernel(method = "rbf", l = 1.25)
  kern_effect_lnr <- 
    parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
  kern_effect_rbf <- 
    parse_kernel_variable("k(x5, x6)", rbf_kern_func, data)
  beta_true <- c(1, .41, 2.37)
  alpha_lnr_true <- rnorm(n)
  alpha_rbf_true <- rnorm(n)
  
  kern_term_lnr <- kern_effect_lnr %*% alpha_lnr_true
  kern_term_rbf <- kern_effect_rbf %*% alpha_rbf_true
  
  data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
    kern_term_lnr + kern_term_rbf
  
  formula_test <- y ~ k(x1):k(x4, x6)
  pvalue <- cvek(formula, kern_func_list = kern_func_list, 
                 data = data, formula_test = formula_test, 
                 alt_kernel_type = "ensemble", B = 200)$pvalue

  expect_lte(pvalue, 1)
  expect_gte(pvalue, 0)
})

test_that("warning message from tuning", {
  kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                         l = c(.5, 1, 1.5), p = 1:3, stringsAsFactors = FALSE)
  # define kernel library
  kern_func_list <- define_library(kern_par)
  n <- 50
  d <- 6
  formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6)
  data <- as.data.frame(matrix(
    rnorm(n * d),
    ncol = d,
    dimnames = list(NULL, paste0("x", 1:d))
  ))
  lnr_kern_func <- generate_kernel(method = "linear")
  rbf_kern_func <- generate_kernel(method = "rbf", l = 1.25)
  kern_effect_lnr <- 
    parse_kernel_variable("k(x3, x4)", lnr_kern_func, data)
  kern_effect_rbf <- 
    parse_kernel_variable("k(x5, x6)", rbf_kern_func, data)
  beta_true <- c(1, .41, 2.37)
  alpha_lnr_true <- rnorm(n)
  alpha_rbf_true <- rnorm(n)
  
  kern_term_lnr <- kern_effect_lnr %*% alpha_lnr_true
  kern_term_rbf <- kern_effect_rbf %*% alpha_rbf_true
  
  data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true + 
    kern_term_lnr + kern_term_rbf
  
  formula_test <- y ~ k(x1):k(x4, x6)
  lambda = rep(.5, 11)
  expect_warning(cvek(formula, kern_func_list = kern_func_list, 
                      data = data, formula_test = formula_test, 
                      lambda = lambda, alt_kernel_type = "ensemble", 
                      B = 200),
                 "the smallest one")
})
