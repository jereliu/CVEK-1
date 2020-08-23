## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
# devtools::install_github("IrisTeng/CVEK-1")
library(CVEK)

## -----------------------------------------------------------------------------
set.seed(0726)
n <- 150 # including training and test
d <- 4
int_effect <- 0.2
data <- matrix(rnorm(n * d), ncol = d)
Z1 <- data[, 1:2]
Z2 <- data[, 3:4]

kern <- generate_kernel(method = "rbf", l = 1)
w <- rnorm(n)
w12 <- rnorm(n)
K1 <- kern(Z1, Z1)
K2 <- kern(Z2, Z2)
K1 <- K1 / sum(diag(K1)) # standardize kernel
K2 <- K2 / sum(diag(K2))
h0 <- K1 %*% w + K2 %*% w
h0 <- h0 / sqrt(sum(h0 ^ 2)) # standardize main effect

h1_prime <- (K1 * K2) %*% w12 # interaction effect

# standardize sampled functions to have unit norm, so that 0.2
# represents the interaction strength relative to main effect
Ks <- svd(K1 + K2)
len <- length(Ks$d[Ks$d / sum(Ks$d) > .001])
U0 <- Ks$u[, 1:len]
h1_prime_hat <- fitted(lm(h1_prime ~ U0))
h1 <- h1_prime - h1_prime_hat

h1 <- h1 / sqrt(sum(h1 ^ 2)) # standardize interaction effect
Y <- h0 + int_effect * h1 + rnorm(1) + rnorm(n, 0, 0.01)
data <- as.data.frame(cbind(Y, Z1, Z2))
colnames(data) <- c("y", paste0("z", 1:d))

data_train <- data[1:100, ]
data_test <- data[101:150, ]

## ---- results='asis'----------------------------------------------------------
knitr::kable(head(data_train, 5))

## ---- fig.width=14, fig.height=11---------------------------------------------
knitr::include_graphics("table1.pdf", auto_pdf = TRUE)

## -----------------------------------------------------------------------------
kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                       l = rep(1, 3), p = 1:3, stringsAsFactors = FALSE)
# define kernel library
kern_func_list <- define_library(kern_par)

## -----------------------------------------------------------------------------
formula <- y ~ z1 + z2 + k(z3, z4)

## -----------------------------------------------------------------------------
est_res <- cvek(formula, kern_func_list = kern_func_list, data = data_train)
est_res$lambda
est_res$u_hat

## -----------------------------------------------------------------------------
formula_test <- y ~ k(z1, z2):k(z3, z4)
pvalue <- cvek(formula, kern_func_list = kern_func_list, 
               data = data_train, formula_test = formula_test,
               mode = "loocv", strategy = "stack",
               beta_exp = 1, lambda = exp(seq(-10, 5)),
               test = "boot", alt_kernel_type = "linear",
               B = 200, verbose = FALSE)$pvalue
pvalue

## -----------------------------------------------------------------------------
y_pred <- predict(est_res, data_test[, 2:5])
data_test_pred <- cbind(y_pred, data_test)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(data_test_pred, 5))

## -----------------------------------------------------------------------------
kern_par <- data.frame(method = c("linear", "polynomial", "rbf"), 
                       l = rep(1, 3), p = 1:3, stringsAsFactors = FALSE)
# define kernel library
kern_func_list <- define_library(kern_par)

formula <- Employed ~ Population + k(Armed.Forces, Population)
formula_test <- Employed ~ k(Armed.Forces):k(Population)
fit_longley <- cvek(formula, kern_func_list = kern_func_list, data = longley, 
                    formula_test = formula_test)
fit_longley$pvalue

