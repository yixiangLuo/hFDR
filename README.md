# Demo version of FDR estimator

## Installation

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("yixiangluo/hFDR")
```

## Simple example with Lasso

``` r
set.seed(2024)

# Problem parameters
n <- 300          # number of observations
p <- 100           # number of variables
k <- 10            # number of variables with nonzero coefficients
amplitude <- 3.5   # signal amplitude (for noise level = 1)

# Generate the variables from a multivariate normal distribution
mu <- rep(0,p)
rho <- 0.25
Sigma <- toeplitz(rho^(0:(p-1)))
X <- matrix(rnorm(n*p),n) %*% chol(Sigma)

# Generate the response from a linear model
nonzero <- sample(p, k)
beta <- amplitude * (1:p %in% nonzero) / sqrt(n)
y.sample <- function(X) X %*% beta + rnorm(n)
y <- y.sample(X)

# generate lamabda sequence for lasso
gene_lambda <- function(X, y){
  n <- NROW(X)
  p <- NCOL(X)
  
  lambda_max <- max(abs(t(X) %*% y)) / n
  sigma <- sqrt(sum(lm(y ~ X)$residuals^2) / (n-p-1))
  lambda_min <- sigma / n
  n_lambda <- 100
  lambda_seq <- lambda_max * (lambda_min/lambda_max)^((0:n_lambda)/n_lambda)
  
  return(lambda_seq)
}
lambda <- gene_lambda(X, y)
```

``` r
library(hFDR)
library(glmnet)
lasso.cv <- cv.glmnet(X, y, lambda = lambda, alpha = 1, nfolds = 10,
                      intercept = T, standardize = T, family = "gaussian")
hFDR.obj <- hFDR(X, y, model = "gausslinear", select = "lasso",
                 lambda = lasso.cv$lambda, nlambda = 40, n_cores = 14)

plot(hFDR.obj, lasso.cv)
```
![hFDR example](readme/hFDR_example.png)
