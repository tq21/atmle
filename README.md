
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`atmle`: Adaptive Targeted Minimum Loss-Based Estimation

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

> This package uses adaptive targeted minimum loss-based estimation to
> estimate the average treatment effect from combined randomized trial
> and real-world data.

**Authors:** [Sky Qiu](https://github.com/tq21), [Mark van der
Laan](https://vanderlaan-lab.org/), [Lars van der
Laan](https://larsvanderlaan.github.io/)

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tq21/atmle/issues).

------------------------------------------------------------------------

## Example

``` r
library(atmle)
library(EScvtmle)
library(sl3)
library(data.table)
library(glmnet)
options(sl3.verbose = TRUE)
`%+%` <- function(a, b) paste0(a, b)
set.seed(82379)

n_rct <- 200 # RCT sample size
n_rwd <- 1000 # RWD sample size
n <- n_rct + n_rwd
g_rct <- 0.67 # RCT treatment probability
true_ate <- 1.5 # true ATE

# error
UY <- rnorm(n, 0, 0.5)
U_bias <- rnorm(n, 0, 0.1)

# baseline covariates
W1 <- rnorm(n, 0, 1)
W2 <- rnorm(n, 0, 1)
W3 <- rnorm(n, 0, 1)
W4 <- rnorm(n, 0, 1)

# study indicator, S=1 for RCT, S=0 for RWD
S <- c(rep(1, n_rct), rep(0, n_rwd))

# treatments (external data has both treated and controls)
A_rct <- rbinom(n_rct, 1, g_rct)
A_rwd <- rbinom(n_rwd, 1, plogis(0.5*W1-0.9*W2))
A <- c(A_rct, A_rwd)

# outcome
Y <- 2.1+0.8*W1+2.5*W2-3.1*W3+0.9*W4+true_ate*A+UY+
  as.numeric(S == 0)*(A*(2.9*W1+2.3*W2+U_bias))

# data frames combining RCT and RWD
data <- data.frame(S = S,
                   W1 = W1,
                   W2 = W2,
                   W3 = W3,
                   W4 = W4,
                   A = A,
                   Y = Y)

# run A-TMLE
res <- atmle(data = data,
             S_node = 1,
             W_node = c(2, 3, 4, 5),
             A_node = 6,
             Y_node = 7,
             controls_only = FALSE,
             atmle_pooled = TRUE,
             var_method = "ic",
             nuisance_method = "sl3",
             working_model = "glmnet",
             g_rct = g_rct,
             verbose = FALSE)
print("A-TMLE ATE estimate: " %+% round(res$est, 2) %+% 
        " (" %+% round(res$lower, 2) %+% ", " %+% round(res$upper, 2) %+% ")")
#> [1] "A-TMLE ATE estimate: 1.55 (1.47, 1.64)"

# compared to ES-CVTMLE
escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3", "W4"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = g_rct,
                          family = "gaussian",
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)
print("ES-CVTMLE ATE estimate: " %+% round(escvtmle_res$ATE$b2v, 2) %+% 
        " (" %+% round(escvtmle_res$CI$b2v[1], 2) %+% 
        ", " %+% round(escvtmle_res$CI$b2v[2], 2) %+% ")")
#> [1] "ES-CVTMLE ATE estimate: 1.57 (1.41, 1.72)"
print("True ATE: " %+% true_ate)
#> [1] "True ATE: 1.5"
```

------------------------------------------------------------------------

## License

© 2023 [Sky Qiu](https://github.com/tq21), [Mark van der
Laan](https://vanderlaan-lab.org/), [Lars van der
Laan](https://larsvanderlaan.github.io/)

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

------------------------------------------------------------------------