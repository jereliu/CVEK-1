## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE----------------------------------------------------
devtools::install_github("IrisTeng/CVEK-1")
library(CVEK)

## ------------------------------------------------------------------------
estimation

## ------------------------------------------------------------------------
set.seed(0726)
n <- 150 # including training and test
d <- 6
data <- as.data.frame(matrix(
  rnorm(n * d),
  ncol = d,
  dimnames = list(NULL, paste0("x", 1:d))
))
beta_true <- c(1, .41, 2.37)
data$y <- as.matrix(cbind(1, data[, c("x1", "x2")])) %*% beta_true

data_train <- data[1:100, ]
data_test <- data[101:150, ]

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(data_train, 10))

## ------------------------------------------------------------------------
kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                         l = c(.5, 1, 1.5), p = 1:3, stringsAsFactors = FALSE)
  # define kernel library
kern_func_list <- define_library(kern_par)

## ------------------------------------------------------------------------
formula <- y ~ x1 + x2 + k(x3, x4) + k(x5, x6)

## ------------------------------------------------------------------------
est_res <- cvek(formula, kern_func_list = kern_func_list, data = data_train)
# est_res$lambda
# est_res$beta
# est_res$u_hat

## ------------------------------------------------------------------------
cvek_test

## ------------------------------------------------------------------------
formula_test <- y ~ k(x1):k(x4, x5)
pvalue <- cvek(formula, kern_func_list = kern_func_list, 
               data = data_train, formula_test = formula_test,
               mode = "loocv", strategy = "stack",
               beta_exp = 1, lambda = exp(seq(-10, 5)),
               test = "boot", alt_kernel_type = "ensemble",
               B = 200, verbose = FALSE)$pvalue

## ------------------------------------------------------------------------
predict.cvek

## ------------------------------------------------------------------------
y_pred <- predict(est_res, data_test[, 1:6])
data_test_pred <- cbind(data_test, y_pred)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(data_test_pred, 10))

