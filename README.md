
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`atmle`: Adaptive Targeted Minimum Loss-Based Estimation

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

> This package uses adaptive targeted minimum loss-based estimation to
> estimate the average treatment effect from combined randomized trial
> and real-world data.

**Authors:** [Sky Qiu](https://github.com/tq21), [Lars van der
Laan](https://larsvanderlaan.github.io/) [Mark van der
Laan](https://vanderlaan-lab.org/),

------------------------------------------------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/tq21/atmle/issues).

------------------------------------------------------------------------

## Example

``` r
library(atmle)
library(EScvtmle)
options(sl3.verbose = TRUE)
`%+%` <- function(a, b) paste0(a, b)
set.seed(82379)

# simulate data
n <- 2000
S <- rbinom(n, 1, 0.2)
W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
A <- numeric(n)
g_rct <- 0.67
A[S == 1] <- rbinom(sum(S), 1, g_rct)
A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
UY <- rnorm(n, 0, 1)
U_bias <- rnorm(n, 0, 0.5)
Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
data <- data.frame(S, W1, W2, A, Y)
true_ate <- 1.5

# run A-TMLE
res <- atmle(data,
             S_node = c(1),
             W_node = c(2, 3),
             A_node = 4,
             Y_node = 5,
             controls_only = FALSE,
             family = "gaussian",
             atmle_pooled = TRUE,
             g_rct = 0.67,
             verbose = FALSE)
print("A-TMLE ATE estimate: " %+% round(res$est, 2) %+% 
        " (" %+% round(res$lower, 2) %+% ", " %+% round(res$upper, 2) %+% ")")
#> [1] "A-TMLE ATE estimate: 1.53 (1.43, 1.62)"
print("True ATE: " %+% true_ate)
#> [1] "True ATE: 1.5"

# compared to ES-CVTMLE
escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2"),
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
#> [1] "ES-CVTMLE ATE estimate: 1.53 (1.33, 1.73)"
print("True ATE: " %+% true_ate)
#> [1] "True ATE: 1.5"
```

------------------------------------------------------------------------

## License

© 2023 [Sky Qiu](https://github.com/tq21), [Lars van der
Laan](https://larsvanderlaan.github.io/) [Mark van der
Laan](https://vanderlaan-lab.org/),

The contents of this repository are distributed under the GPL-3 license.
See file `LICENSE` for details.

------------------------------------------------------------------------
