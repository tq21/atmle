test_that("learn_theta_tilde works when external has both treated and controls,
          and outcome is continuous", {
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

  theta_tilde <- learn_theta_tilde(W, Y, "glm", 5, "gaussian", c(0, 0.5)) # glm
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 0.5))

  theta_tilde <- learn_theta_tilde(W, Y, "glmnet", 5, "gaussian", c(-4, 4)) # glmnet
  expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))

  theta_tilde <- learn_theta_tilde(W, Y, "sl3", 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta_tilde <- learn_theta_tilde(W, Y, lrnrs, 5, "gaussian", c(-4, 4)) # sl3 learners
  expect_true(all(theta_tilde >= -4 & theta_tilde <= 4))

  expect_error(learn_theta_tilde(W, Y, "abc", 5, "gaussian", c(-4, 4)))
})

test_that("learn_theta_tilde works when external has both treated and controls,
          and outcome is binary", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
  Y <- rbinom(n, 1, plogis(-0.5-0.8*W1-1.1*W2+1.5*A+(1-S)*(0.9+2.6*W1)))

  theta_tilde <- learn_theta_tilde(W, Y, "glm", 5, "binomial", c(0.2, 0.8)) # glm
  expect_true(all(theta_tilde >= 0.2 & theta_tilde <= 0.8))

  theta_tilde <- learn_theta_tilde(W, Y, "glmnet", 5, "binomial", NULL) # glmnet
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))

  theta_tilde <- learn_theta_tilde(W, Y, "sl3", 5, "binomial", NULL) # sl3
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  theta_tilde <- learn_theta_tilde(W, Y, lrnrs, 5, "binomial", NULL) # sl3 learners
  expect_true(all(theta_tilde >= 0 & theta_tilde <= 1))

  expect_error(learn_theta_tilde(W, Y, "abc", 5, "binomial", NULL))
})
