[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/robStepSplitReg)](https://cran.r-project.org/package=robStepSplitReg)
[![CRAN Data](https://www.r-pkg.org/badges/last-release/robStepSplitReg)](https://cran.r-project.org/package=robStepSplitReg) 
[![Downloads](http://cranlogs.r-pkg.org/badges/robStepSplitReg)](https://cran.r-project.org/package=robStepSplitReg)

robStepSplitReg
================

This package provides functions for performing robust stepwise split regularized regression.

---------------------------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=robStepSplitReg).

```{r installation, eval = FALSE}
install.packages("robStepSplitReg", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/robStepSplitReg)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/robStepSplitReg")
```

### Usage

``` r
# Required library
library(mvnfast)

# Simulation parameters
n <- 50
p <- 500
rho <- 0.8
p.active <- 100
snr <- 3
contamination.prop <- 0.2

# Setting the seed
set.seed(0)

# Simulation of beta vector
true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))

# Simulation of uncontaminated data 
sigma.mat <- matrix(0, nrow = p, ncol = p)
sigma.mat[1:p.active, 1:p.active] <- rho
diag(sigma.mat) <- 1
x <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))
y <- x %*% true.beta + rnorm(n, 0, sigma)

# Contamination of data 
contamination_indices <- 1:floor(n*contamination.prop)
k_lev <- 2
k_slo <- 100
x_train <- x
y_train <- y
beta_cont <- true.beta
beta_cont[true.beta!=0] <- beta_cont[true.beta!=0]*(1 + k_slo)
beta_cont[true.beta==0] <- k_slo*max(abs(true.beta))
for(cont_id in contamination_indices){
  
  a <- runif(p, min = -1, max = 1)
  a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
  x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
    k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
  y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
}

# Ensemble models
ensemble_fit <- robStepSplitReg(x_train, y_train,
                                n_models = 5,
                                model_saturation = c("fixed", "p-value")[1],
                                alpha = 0.05, model_size = 25,
                                robust = TRUE,
                                compute_coef = TRUE,
                                pense_alpha = 1/4, pense_cv_k = 5, pense_cv_repl = 1,
                                cl = NULL)

# Ensemble coefficients
ensemble_coefs <- coef(ensemble_fit, group_index = 1:ensemble_fit$n_models)
sens_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/p.active
spec_ensemble <- sum(which((ensemble_coefs[-1]!=0)) <= p.active)/sum(ensemble_coefs[-1]!=0)

# Simulation of test data
m <- 2e3
x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test %*% true.beta + rnorm(m, 0, sigma)

# Prediction of test samples
ensemble_preds <- predict(ensemble_fit, newx = x_test, 
                          group_index = 1:ensemble_fit$n_models,
                          dynamic = FALSE)
mspe_ensemble <- mean((y_test - ensemble_preds)^2)/sigma^2
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
