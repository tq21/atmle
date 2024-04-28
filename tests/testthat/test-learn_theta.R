test_that("learn_theta works when external has both treated and controls, and
          outcome is continuous", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1[S == 0]-0.9*W2[S == 0]))
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # glm
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = FALSE,
                       method = "glm",
                       folds = folds,
                       family = "gaussian",
                       theta_bounds = c(0, 0.5))
  expect_true(all(theta >= 0 & theta <= 0.5))
  expect_equal(length(theta), n)

  # glmnet
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = FALSE,
                       method = "glmnet",
                       folds = folds,
                       family = "gaussian",
                       theta_bounds = c(-4, 4))
  expect_true(all(theta >= -4 & theta <= 4))
  expect_equal(length(theta), n)

  # # sl3 with default learners
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = FALSE,
  #                      method = "sl3",
  #                      folds = folds,
  #                      family = "gaussian",
  #                      theta_bounds = c(-4, 4)) # sl3
  # expect_true(all(theta >= -4 & theta <= 4))
  # expect_equal(length(theta), n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = FALSE,
  #                      method = lrnrs,
  #                      folds = folds,
  #                      family = "gaussian",
  #                      theta_bounds = c(-4, 4))
  # expect_true(all(theta >= -4 & theta <= 4))
  # expect_equal(length(theta), n)

  expect_error(learn_theta(W = W,
                           A = A,
                           Y = Y,
                           delta = delta,
                           controls_only = FALSE,
                           method = "abc",
                           folds = folds,
                           family = "gaussian",
                           theta_bounds = c(-4, 4)))
})

test_that("learn_theta works when external has both treated and controls, and
          outcome is binary", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1[S == 0]-0.9*W2[S == 0]))
  Y <- rbinom(n, 1, plogis(-0.5-0.8*W1-1.1*W2+1.5*A+(1-S)*(0.9+2.6*W1)))
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # glm
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = FALSE,
                       method = "glm",
                       folds = folds,
                       family = "binomial",
                       theta_bounds = c(0.2, 0.8))
  expect_true(all(theta >= 0.2 & theta <= 0.8))
  expect_equal(length(theta), n)

  # glmnet
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = FALSE,
                       method = "glm",
                       folds = folds,
                       family = "binomial",
                       theta_bounds = NULL)
  expect_true(all(theta >= 0 & theta <= 1))
  expect_equal(length(theta), n)

  # # sl3 with default learners
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = FALSE,
  #                      method = "sl3",
  #                      folds = folds,
  #                      family = "binomial",
  #                      theta_bounds = NULL)
  # expect_true(all(theta >= 0 & theta <= 1))
  # expect_equal(length(theta), n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = FALSE,
  #                      method = lrnrs,
  #                      folds = folds,
  #                      family = "binomial",
  #                      theta_bounds = NULL)
  # expect_true(all(theta >= 0 & theta <= 1))
  # expect_equal(length(theta), n)

  expect_error(learn_theta(W = W,
                           A = A,
                           Y = Y,
                           delta = delta,
                           controls_only = FALSE,
                           method = "abc",
                           folds = folds,
                           family = "binomial",
                           theta_bounds = NULL))
})

test_that("learn_theta works when external has controls only, and outcome
          is continuous", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0
  UY <- rnorm(n, 0, 1)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = TRUE,
                       method = "glm",
                       folds = folds,
                       family = "gaussian",
                       theta_bounds = c(0, 0.5)) # glm
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 0.5))
  expect_length(theta, n)

  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = TRUE,
                       method = "glmnet",
                       folds = folds,
                       family = "gaussian",
                       theta_bounds = c(-4, 4))
  expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  expect_length(theta, n)

  # # sl3 with default learners
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = TRUE,
  #                      method = "sl3",
  #                      folds = folds,
  #                      family = "gaussian",
  #                      theta_bounds = c(-4, 4))
  # expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  # expect_length(theta, n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = TRUE,
  #                      method = lrnrs,
  #                      folds = folds,
  #                      family = "gaussian",
  #                      theta_bounds = c(-4, 4))
  # expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  # expect_length(theta, n)

  expect_error(learn_theta(W = W,
                           A = A,
                           Y = Y,
                           delta = delta,
                           controls_only = TRUE,
                           method = "abc",
                           folds = folds,
                           family = "gaussian",
                           theta_bounds = c(-4, 4)))
})

test_that("learn_theta works when external has controls only, and outcome
          is binary", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0
  Y <- rbinom(n, 1, plogis(-0.5-0.8*W1-1.1*W2+1.5*A+(1-S)*(0.9+2.6*W1)))
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(n = n, V = 5,
                        strata_ids = as.integer(factor(cv_strata)))
  })

  # glm
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = TRUE,
                       method = "glm",
                       folds = folds,
                       family = "binomial",
                       theta_bounds = c(0.2, 0.8)) # glm
  expect_true(all(theta[A == 0] >= 0.2 & theta[A == 0] <= 0.8))
  expect_length(theta, n)

  # glmnet
  theta <- learn_theta(W = W,
                       A = A,
                       Y = Y,
                       delta = delta,
                       controls_only = TRUE,
                       method = "glmnet",
                       folds = folds,
                       family = "binomial",
                       theta_bounds = NULL)
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  expect_length(theta, n)

  # # sl3 with default learners
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = TRUE,
  #                      method = "sl3",
  #                      folds = folds,
  #                      family = "binomial",
  #                      theta_bounds = NULL)
  # expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  # expect_length(theta, n)
#
  # # sl3 with custom learners
  # lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  # theta <- learn_theta(W = W,
  #                      A = A,
  #                      Y = Y,
  #                      delta = delta,
  #                      controls_only = TRUE,
  #                      method = lrnrs,
  #                      folds = folds,
  #                      family = "binomial",
  #                      theta_bounds = NULL) # sl3
  # expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  # expect_length(theta, n)

  expect_error(learn_theta(W = W,
                           A = A,
                           Y = Y,
                           delta = delta,
                           controls_only = TRUE,
                           method = "abc",
                           folds = folds,
                           family = "binomial",
                           theta_bounds = NULL))
})
