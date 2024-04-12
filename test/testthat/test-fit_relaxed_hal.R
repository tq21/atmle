test_that("basic functionality test", {
  # generate data
  set.seed(123)
  n <- 100; p <- 2
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  Y <- X[, 1] + X[, 2] + rnorm(n)

  # fit
  fit <- fit_relaxed_hal(X = X,
                         Y = Y,
                         X_unpenalized = NULL,
                         X_weak_penalized = NULL,
                         X_weak_penalized_level = 0,
                         family = "gaussian",
                         weights = NULL,
                         relaxed = TRUE,
                         v_folds = 5)

  expect_equal(dim(fit$x_basis)[1], length(fit$pred))
  expect_equal(dim(fit$x_basis)[2], length(fit$beta),
               length(fit$hal_basis_list))
})
