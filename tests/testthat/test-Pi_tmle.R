library(atmle)
source("tests/testthat/utils.R")

test_Pi_tmle <- function(test_args) {
  set.seed(123)
  data <- sim_data(test_args$n,
                   test_args$controls_only,
                   test_args$family,
                   test_args$prop_miss)
  S <- data$S; W <- data[, c("W1", "W2")]; A <- data$A; Y <- data$Y
  delta <- as.integer(!is.na(Y))
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = test_args$n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = test_args$controls_only,
                       method = test_args$theta_method,
                       folds = folds,
                       family = test_args$family,
                       theta_bounds = test_args$theta_bounds,
                       cross_fit_nuisance = test_args$cross_fit_nuisance)

  g <- learn_g(W = W,
               A = A,
               method = test_args$g_method,
               folds = folds,
               g_bounds = test_args$g_bounds)

  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = test_args$controls_only,
                 method = test_args$Pi_method,
                 folds = folds,
                 Pi_bounds = test_args$Pi_bounds,
                 cross_fit_nuisance = test_args$cross_fit_nuisance)

  tau <- learn_tau(S = S, W = W, A = A, Y = Y, Pi = Pi, theta = theta, g = g,
                   delta = delta,
                   controls_only = test_args$controls_only,
                   method = test_args$bias_working_model,
                   v_folds = length(folds),
                   max_degree = test_args$max_degree,
                   min_working_model = test_args$min_working_model,
                   target_gwt = test_args$target_gwt,
                   Pi_bounds = test_args$Pi_bounds,
                   enumerate_basis_args = test_args$enumerate_basis_args,
                   fit_hal_args = test_args$fit_hal_args,
                   weights = test_args$weights,
                   bias_working_model_formula = test_args$bias_working_model_formula)

  res <- Pi_tmle(S = S,
                 W = W,
                 A = A,
                 g = g,
                 tau = tau,
                 Pi = Pi,
                 controls_only = test_args$controls_only,
                 target_gwt = test_args$target_gwt,
                 Pi_bounds = test_args$Pi_bounds)

  return(res)
}

test_that("Pi_tmle basic functionalities", {
  # tabulate test scenarios
  test_args <- list.expand(n = list(500),
                           controls_only = list(TRUE, FALSE),
                           family = list("gaussian", "binomial"),
                           prop_miss = list(0, 0.1),
                           theta_method = list("glm"),
                           theta_bounds = list(c(-Inf, Inf)),
                           g_method = list("glm"),
                           g_bounds = list(c(0.01, 0.99)),
                           Pi_method = list("glm"),
                           Pi_bounds = list(c(0.1, 0.9)),
                           bias_working_model = list("glmnet"),
                           max_degree = list(1),
                           min_working_model = list(FALSE),
                           target_gwt = list(TRUE, FALSE),
                           enumerate_basis_args = list(list()),
                           fit_hal_args = list(list()),
                           weights = list(rep(1, 500)),
                           bias_working_model_formula = list(NULL),
                           cross_fit_nuisance = list(TRUE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_Pi_tmle(cur_args)

    # checks
    expect_length(res$pred, cur_args$n)
    expect_true(all(res$A0 >= cur_args$Pi_bounds[1] & res$A0 <= cur_args$Pi_bounds[2]))
    if (cur_args$controls_only) {
      expect_true(all(res$A1 == 1))
    } else {
      expect_true(all(res$A1 >= cur_args$Pi_bounds[1] & res$A1 <= cur_args$Pi_bounds[2]))
    }
  }

  # TODO: missing data case
})
