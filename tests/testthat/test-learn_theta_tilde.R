test_that("learn_theta_tilde works when outcome is continuous", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n - sum(S), 1, plogis(1.2 * W1[S == 0] - 0.9 * W2[S == 0]))
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5 - 0.8 * W1 - 1.1 * W2 + 1.5 * A + UY + (1 - S) * (0.9 + 2.6 * W1) + (1 - S) * U_bias
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # glm
  theta_tilde <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "glm",
    folds = folds,
    family = "gaussian",
    theta_bounds = c(0, 0.5)
  )
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 0.5))
  expect_length(theta_tilde, n)

  # glmnet
  theta_tilde <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "glmnet",
    folds = folds,
    family = "gaussian",
    theta_bounds = c(-4, 4)
  )
  expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))
  expect_length(theta_tilde, n)

  # # sl3 with default learners
  # theta_tilde <- learn_theta_tilde(W = W,
  #                                  Y = Y,
  #                                  delta = delta,
  #                                  method = "sl3",
  #                                  folds = folds,
  #                                  family = "gaussian",
  #                                  theta_bounds = c(-4, 4))
  # expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))
  # expect_length(theta_tilde, n)
  #
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta_tilde <- learn_theta_tilde(W = W,
  #                                  Y = Y,
  #                                  delta = delta,
  #                                  method = lrnrs,
  #                                  folds = folds,
  #                                  family = "gaussian",
  #                                  theta_bounds = c(-4, 4))
  # expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))
  # expect_length(theta_tilde, n)

  expect_error(learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "abc",
    folds = folds,
    family = "gaussian",
    theta_bounds = c(-4, 4)
  ))
})

test_that("learn_theta_tilde works when outcome is binary", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n - sum(S), 1, plogis(1.2 * W1[S == 0] - 0.9 * W2[S == 0]))
  Y <- rbinom(n, 1, plogis(-0.5 - 0.8 * W1 - 1.1 * W2 + 1.5 * A + (1 - S) * (0.9 + 2.6 * W1)))
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # glm
  theta_tilde <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "glm",
    folds = folds,
    family = "binomial",
    theta_bounds = c(0.2, 0.5)
  )
  expect_true(all(theta_tilde >= 0.2 & theta_tilde <= 0.5))
  expect_length(theta_tilde, n)

  # glm
  theta_tilde <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "glmnet",
    folds = folds,
    family = "binomial",
    theta_bounds = NULL
  )
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))
  expect_length(theta_tilde, n)

  # sl3 with default learners
  # theta_tilde <- learn_theta_tilde(W = W,
  #                                  Y = Y,
  #                                  delta = delta,
  #                                  method = "sl3",
  #                                  folds = folds,
  #                                  family = "binomial",
  #                                  theta_bounds = NULL)
  # expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))
  # expect_length(theta_tilde, n)
  #
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta_tilde <- learn_theta_tilde(W = W,
  #                                  Y = Y,
  #                                  delta = delta,
  #                                  method = lrnrs,
  #                                  folds = folds,
  #                                  family = "binomial",
  #                                  theta_bounds = NULL)
  # expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))
  # expect_length(theta_tilde, n)

  expect_error(learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = "abc",
    folds = folds,
    family = "binomial",
    theta_bounds = NULL
  ))
})
