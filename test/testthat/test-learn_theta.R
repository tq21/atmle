test_that("learn_theta works when external has both treated and controls, and
          outcome is continuous", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias

  theta <- learn_theta(W, A, Y, FALSE, "glm", 5, "gaussian", c(0, 0.5)) # glm
  expect_true(all(theta >= 0 & theta <= 0.5))

  theta <- learn_theta(W, A, Y, FALSE, "glmnet", 5, "gaussian", c(-4, 4)) # glmnet
  expect_true(all(theta >= -4 & theta <= 4))

  theta <- learn_theta(W, A, Y, FALSE, "sl3", 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(theta >= -4 & theta <= 4))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta <- learn_theta(W, A, Y, FALSE, lrnrs, 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(theta >= -4 & theta <= 4))

  expect_error(learn_theta(W, A, Y, FALSE, "abc", 5, "gaussian", c(-4, 4)))
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
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
  Y <- rbinom(n, 1, plogis(-0.5-0.8*W1-1.1*W2+1.5*A+(1-S)*(0.9+2.6*W1)))

  theta <- learn_theta(W, A, Y, FALSE, "glm", 5, "binomial", c(0.2, 0.8)) # glm
  expect_true(all(theta >= 0.2 & theta <= 0.8))

  theta <- learn_theta(W, A, Y, FALSE, "glmnet", 5, "binomial", NULL) # glmnet
  expect_true(all(theta >= 0 & theta <= 1))

  theta <- learn_theta(W, A, Y, FALSE, "sl3", 5, "binomial", NULL) # sl3
  expect_true(all(theta >= 0 & theta <= 1))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta <- learn_theta(W, A, Y, FALSE, lrnrs, 5, "binomial", NULL) # sl3
  expect_true(all(theta >= 0 & theta <= 1))

  expect_error(learn_theta(W, A, Y, FALSE, "abc", 5, "binomial", NULL))
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
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias

  theta <- learn_theta(W, A, Y, TRUE, "glm", 5, "gaussian", c(0, 0.5)) # glm
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 0.5))
  expect_true(all(theta[A == 1] == 0))

  theta <- learn_theta(W, A, Y, TRUE, "glmnet", 5, "gaussian", c(-4, 4)) # glmnet
  expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  expect_true(all(theta[A == 1] == 0))

  theta <- learn_theta(W, A, Y, TRUE, "sl3", 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  expect_true(all(theta[A == 1] == 0))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta <- learn_theta(W, A, Y, TRUE, lrnrs, 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(theta[A == 0] >= -4 & theta[A == 0] <= 4))
  expect_true(all(theta[A == 1] == 0))

  expect_error(learn_theta(W, A, Y, FALSE, "abc", 5, "gaussian", c(-4, 4)))
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

  theta <- learn_theta(W, A, Y, TRUE, "glm", 5, "binomial", c(0.2, 0.8)) # glm
  expect_true(all(theta[A == 0] >= 0.2 & theta[A == 0] <= 0.8))
  expect_true(all(theta[A == 1] == 0))

  theta <- learn_theta(W, A, Y, TRUE, "glmnet", 5, "binomial", NULL) # glmnet
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  expect_true(all(theta[A == 1] == 0))

  theta <- learn_theta(W, A, Y, TRUE, "sl3", 5, "binomial", NULL) # sl3
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  expect_true(all(theta[A == 1] == 0))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta <- learn_theta(W, A, Y, TRUE, lrnrs, 5, "binomial", NULL) # sl3
  expect_true(all(theta[A == 0] >= 0 & theta[A == 0] <= 1))
  expect_true(all(theta[A == 1] == 0))

  expect_error(learn_theta(W, A, Y, FALSE, "abc", 5, "binomial", NULL))
})
