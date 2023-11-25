test_that("Pi_tmle works when external has both treated and controls", {
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

  g <- learn_g(S, W, A, 0.67, FALSE, "glm", 5, c(0.01, 0.99)) # glm
  Pi <- learn_Pi(S, W, A, FALSE, "glm", 5, c(0.01, 0.99))
  theta <- learn_theta(W, A, Y, FALSE, "glm", 5, "gaussian", NULL)
  tau <- learn_tau(S, W, A, Y, Pi, theta, FALSE, "glmnet", 5, 0)
  Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, FALSE, TRUE, c(0.01, 0.99)) # target_gwt = TRUE
  expect_equal(n, length(Pi_star$pred))
  expect_equal(n, length(Pi_star$A0))
  expect_equal(n, length(Pi_star$A1))
  expect_true(all(Pi_star$pred >= 0.01 & Pi_star$pred <= 0.99))
  expect_true(all(Pi_star$A0 >= 0.01 & Pi_star$A0 <= 0.99))
  expect_true(all(Pi_star$A1 >= 0.01 & Pi_star$A1 <= 0.99))

  Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, FALSE, FALSE, c(0.01, 0.99)) # target_gwt = FALSE
  expect_equal(n, length(Pi_star$pred))
  expect_equal(n, length(Pi_star$A0))
  expect_equal(n, length(Pi_star$A1))
  expect_true(all(Pi_star$pred >= 0.01 & Pi_star$pred <= 0.99))
  expect_true(all(Pi_star$A0 >= 0.01 & Pi_star$A0 <= 0.99))
  expect_true(all(Pi_star$A1 >= 0.01 & Pi_star$A1 <= 0.99))
})

test_that("Pi_tmle works when external has controls only", {
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

  g <- learn_g(S, W, A, 0.67, TRUE, "glm", 5, c(0.01, 0.99)) # glm
  Pi <- learn_Pi(S, W, A, TRUE, "glm", 5, c(0.01, 0.99))
  theta <- learn_theta(W, A, Y, TRUE, "glm", 5, "gaussian", NULL)
  tau <- learn_tau(S, W, A, Y, Pi, theta, TRUE, "glmnet", 5, 0)
  Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, TRUE, TRUE, c(0.01, 0.99)) # target_gwt = TRUE
  expect_equal(n, length(Pi_star$pred))
  expect_equal(n, length(Pi_star$A0))
  expect_true(all(Pi_star$A1 == 1))
  expect_true(all(Pi_star$pred[A == 1] == 1))
  expect_true(all(Pi_star$pred[A == 0] >= 0.01 & Pi_star$pred[A == 0] <= 0.99))
  expect_true(all(Pi_star$A0 >= 0.01 & Pi_star$A0 <= 0.99))

  Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, TRUE, FALSE, c(0.01, 0.99)) # target_gwt = FALSE
  expect_equal(n, length(Pi_star$pred))
  expect_equal(n, length(Pi_star$A0))
  expect_true(all(Pi_star$A1 == 1))
  expect_true(all(Pi_star$pred[A == 1] == 1))
  expect_true(all(Pi_star$pred[A == 0] >= 0.01 & Pi_star$pred[A == 0] <= 0.99))
  expect_true(all(Pi_star$A0 >= 0.01 & Pi_star$A0 <= 0.99))
})
