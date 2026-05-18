
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `atmle`: Adaptive TMLE for RCT + RWD Data Fusion

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

This package exposes the `atmle_ate_fusion` R6 estimator for adaptive
targeted minimum loss-based estimation of average treatment effects using
combined randomized trial and real-world data.

## Example

``` r
library(atmle)
library(sl3)

set.seed(147)
n <- 160
W1 <- rnorm(n)
W2 <- rnorm(n)
S <- rbinom(n, 1, plogis(-0.3 + 0.25 * W1))
A <- rbinom(n, 1, plogis(0.1 * S + 0.2 * W1 - 0.1 * W2))
Y <- 0.5 + 0.4 * W1 + 0.2 * W2 + A + 0.2 * (1 - S) * W2 + rnorm(n)
data <- data.frame(S, W1, W2, A, Y)

method <- list(learners = list(Lrnr_glm$new()))

fit <- atmle_ate_fusion$new(
  data = data,
  S_node = "S",
  W_nodes = c("W1", "W2"),
  A_node = "A",
  Y_node = "Y",
  family = "gaussian",
  n_folds = 3
)

fit$run(
  g_bar_method = method,
  theta_method = method,
  Pi_method = method,
  Q_bar_method = method,
  A_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 5L),
  S_cate_args = list(max_degree = 1L, smoothness_orders = 1L, num_knots = 5L),
  target_method = "tmle",
  target_gwt = FALSE,
  max_iter = 1,
  n_lambda = 1,
  verbose = FALSE
)

fit$results
```

## License

The contents of this repository are distributed under the GPL-3 license.
