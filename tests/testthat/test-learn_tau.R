library(atmle)
source("tests/testthat/utils.R")

test_learn_tau <- function(test_args) {
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

  Pi <- learn_Pi(S = S,
                 W = W,
                 A = A,
                 controls_only = test_args$controls_only,
                 method = test_args$method,
                 folds = folds,
                 Pi_bounds = test_args$Pi_bounds,
                 cross_fit_nuisance = test_args$cross_fit_nuisance)

  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = test_args$controls_only,
                       method = test_args$method,
                       folds = folds,
                       family = test_args$family,
                       theta_bounds = test_args$theta_bounds,
                       cross_fit_nuisance = test_args$cross_fit_nuisance)

  res <- learn_tau(S = S, W = W, A = A, Y = Y, Pi = Pi, theta = theta, g = g,
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

  return(res)
}

test_that("learn_tau basic functionalities with glmnet working model", {
  # tabulate test scenarios
  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  test_args <- list.expand(n = list(500),
                           controls_only = list(TRUE, FALSE),
                           family = list("gaussian", "binomial"),
                           prop_miss = list(0, 0.1),
                           method = list("glm"),
                           bias_working_model = list("glmnet"),
                           max_degree = list(1),
                           min_working_model = list(FALSE),
                           target_gwt = list(TRUE),
                           Pi_bounds = list(c(0.01, 0.99)),
                           theta_bounds = list(c(-Inf, Inf)),
                           g_bounds = list(c(0.01, 0.99)),
                           enumerate_basis_args = list(list()),
                           fit_hal_args = list(list()),
                           weights = list(rep(1, 500)),
                           bias_working_model_formula = list(NULL),
                           cross_fit_nuisance = list(TRUE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_learn_tau(cur_args)

    # checks
    expect_length(res$A1, cur_args$n); expect_length(res$A0, cur_args$n)
    expect_equal(dim(res$x_basis_A0)[1], cur_args$n)
    expect_equal(as.numeric(res$x_basis%*%matrix(res$coefs)), res$pred)
    expect_equal(as.numeric(res$x_basis_A0%*%matrix(res$coefs)), res$A0)
    expect_length(res$coefs, length(res$non_zero))
    if (cur_args$controls_only) {
      expect_true(all(res$A1 == 0))
      expect_equal(dim(res$x_basis)[1], cur_args$n)
      expect_null(res$x_basis_A1)
    } else {
      expect_equal(as.numeric(res$x_basis_A1%*%matrix(res$coefs)), res$A1)
    }
  }

  # TODO: missing data case
})
