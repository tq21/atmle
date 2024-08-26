library(atmle)
source("tests/testthat/utils.R")

test_learn_T <- function(test_args) {
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

  g <- learn_g(W = W,
               A = A,
               method = test_args$g_method,
               folds = folds,
               g_bounds = test_args$g_bounds,
               cross_fit_nuisance = test_args$cross_fit_nuisance)

  theta_tilde <- learn_theta_tilde(W = W,
                                   Y = Y,
                                   delta = delta,
                                   method = test_args$theta_tilde_method,
                                   folds = folds,
                                   family = test_args$family,
                                   theta_bounds = test_args$theta_bounds,
                                   cross_fit_nuisance = test_args$cross_fit_nuisance)

  res <- learn_T(W = W, A = A, Y = Y, g = g, delta = delta, theta_tilde = theta_tilde,
                 method = test_args$pooled_working_model,
                 v_folds = length(folds),
                 weights = test_args$weights,
                 enumerate_basis_args = test_args$enumerate_basis_args,
                 fit_hal_args = test_args$fit_hal_args,
                 pooled_working_model_formula = test_args$pooled_working_model_formula)

  return(res)
}

test_that("learn_T basic functionalities with glmnet working model", {
  # tabulate test scenarios
  test_args <- list.expand(n = list(500),
                           controls_only = FALSE,
                           family = list("gaussian", "binomial"),
                           prop_miss = list(0, 0.1),
                           g_method = list("glm"),
                           g_bounds = list(c(0.01, 0.99)),
                           theta_tilde_method = list("glm"),
                           theta_bounds = list(c(-Inf, Inf)),
                           pooled_working_model = list("glmnet"),
                           weights = list(rep(1, 500)),
                           enumerate_basis_args = list(list()),
                           fit_hal_args = list(list()),
                           pooled_working_model_formula = list(NULL),
                           cross_fit_nuisance = list(TRUE))

  for (i in 1:length(test_args)) {
    cur_args <- test_args[[i]]
    res <- test_learn_T(cur_args)

    # checks
    expect_length(res$pred, cur_args$n)
    expect_equal(dim(res$x_basis)[1], cur_args$n)
    expect_equal(dim(res$x_basis)[2], length(res$non_zero))
    expect_equal(as.numeric(res$x_basis%*%matrix(res$coefs)), res$pred)
  }

  # TODO: missing data case
})
