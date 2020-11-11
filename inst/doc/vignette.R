## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE---------------------------------------------------------
# devtools::install_github("IrisTeng/CVEK-1")
library(CVEK)
library(ggplot2)
library(ggrepel)

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
est_res <- cvek(formula, kern_func_list = kern_func_list, 
                data = data_train)
est_res$lambda
est_res$u_hat

## -----------------------------------------------------------------------------
formula_test <- y ~ k(z1, z2):k(z3, z4)

cvek(formula, kern_func_list = kern_func_list, 
     data = data_train, formula_test = formula_test,
     mode = "loocv", strategy = "stack",
     beta_exp = 1, lambda = exp(seq(-10, 5)),
     test = "asymp", alt_kernel_type = "ensemble",
     verbose = FALSE)$pvalue

cvek(formula, kern_func_list = kern_func_list, 
     data = data_train, formula_test = formula_test,
     mode = "loocv", strategy = "stack",
     beta_exp = 1, lambda = exp(seq(-10, 5)),
     test = "boot", alt_kernel_type = "ensemble",
     B = 200, verbose = FALSE)$pvalue

## -----------------------------------------------------------------------------
y_pred <- predict(est_res, data_test[, 2:5])
data_test_pred <- cbind(y_pred, data_test)

## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(data_test_pred, 5))

## ---- fig.width=14, fig.height=3----------------------------------------------
knitr::include_graphics("table2.pdf", auto_pdf = TRUE)

## -----------------------------------------------------------------------------
kern_par <- data.frame(method = c("linear", "rbf"), 
                       l = rep(1, 2), p = 1:2, stringsAsFactors = FALSE)
# define kernel library
kern_func_list <- define_library(kern_par)

## -----------------------------------------------------------------------------
formula <- medv ~ zn + indus + chas + nox + rm + age + dis + 
  rad + tax + ptratio + black + k(crim) + k(lstat)
formula_test <- medv ~ k(crim):k(lstat)
fit_bos <- cvek(formula, kern_func_list = kern_func_list, data = Boston, 
                formula_test = formula_test, 
                lambda = exp(seq(-3, 5)), test = "asymp")

## -----------------------------------------------------------------------------
fit_bos$pvalue

## ----message=FALSE, fig.width=7, fig.height=5---------------------------------
# first fit the alternative model
formula_alt <- medv ~ zn + indus + chas + nox + rm + age + dis + 
  rad + tax + ptratio + black + k(crim):k(lstat)
fit_bos_alt <- cvek(formula = formula_alt, kern_func_list = kern_func_list, 
                    data = Boston, lambda = exp(seq(-3, 5)))

# mean-center all confounding variables not involved in the interaction 
# so that the predicted values are more easily interpreted
pred_name <- c("zn", "indus", "chas", "nox", "rm", "age", 
               "dis", "rad", "tax", "ptratio", "black")
covar_mean <- apply(Boston, 2, mean)
pred_cov <- covar_mean[pred_name]
pred_cov_df <- t(as.data.frame(pred_cov))
lstat_list <- seq(12.5, 17.5, length.out = 100)
crim_quantiles <- quantile(Boston$crim, probs = c(.05, .25, .5, .75, .95))

# crim is set to its 5% quantile
data_test1 <- data.frame(pred_cov_df, lstat = lstat_list, 
                             crim = crim_quantiles[1])
data_test1_pred <- predict(fit_bos_alt, data_test1)

# crim is set to its 25% quantile
data_test2 <- data.frame(pred_cov_df, lstat = lstat_list, 
                             crim = crim_quantiles[2])
data_test2_pred <- predict(fit_bos_alt, data_test2)

# crim is set to its 50% quantile
data_test3 <- data.frame(pred_cov_df, lstat = lstat_list, 
                             crim = crim_quantiles[3])
data_test3_pred <- predict(fit_bos_alt, data_test3)

# crim is set to its 75% quantile
data_test4 <- data.frame(pred_cov_df, lstat = lstat_list, 
                             crim = crim_quantiles[4])
data_test4_pred <- predict(fit_bos_alt, data_test4)

# crim is set to its 95% quantile
data_test5 <- data.frame(pred_cov_df, lstat = lstat_list, 
                             crim = crim_quantiles[5])
data_test5_pred <- predict(fit_bos_alt, data_test5)

# combine five sets of prediction data together
medv <- rbind(data_test1_pred, data_test2_pred, data_test3_pred, 
              data_test4_pred, data_test5_pred)
data_pred <- data.frame(lstat = rep(lstat_list, 5), medv = medv, 
                        crim = rep(c("5% quantile", "25% quantile", 
                                     "50% quantile", "75% quantile", 
                                     "95% quantile"), each = 100))
data_pred$crim <- factor(data_pred$crim, 
                         levels = c("5% quantile", "25% quantile", 
                                    "50% quantile", "75% quantile", 
                                    "95% quantile"))

data_label <- data_pred[which(data_pred$lstat == 17.5), ]
data_label$value <- c("0.028%", "0.082%", "0.257%", "3.677%", "15.789%")
data_label$value <- factor(data_label$value, levels = 
                             c("0.028%", "0.082%", "0.257%", 
                               "3.677%", "15.789%"))

ggplot(data = data_pred, aes(x = lstat, y = medv, color = crim)) + 
  geom_point(size = 0.1) + 
  geom_text_repel(aes(label = value), data = data_label, 
                  color = "black", size = 3.6) + 
    scale_colour_manual(values = c("firebrick1", "chocolate2", 
                                   "darkolivegreen3", "skyblue2", 
                                   "purple2")) + 
  geom_line() + theme_set(theme_bw()) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12)) + 
    labs(x = "percentage of lower status", 
         y = "median value of owner-occupied homes ($1000)", 
         col = "per capita crime rate")

