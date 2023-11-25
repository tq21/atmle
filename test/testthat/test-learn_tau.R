test_that("learn_tau works when external has both treated and controls", {
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

  Pi <- learn_Pi(S, W, A, FALSE, "glm", 5, c(0.01, 0.99))
  theta <- learn_theta(W, A, Y, FALSE, "glm", 5, "gaussian", NULL)
  tau <- learn_tau(S, W, A, Y, Pi, theta, FALSE, "glmnet", 5, 0)

  expect_equal(length(tau$A0), n)
  expect_equal(length(tau$A1), n)
  expect_equal(length(tau$pred), n)
  expect_equal(dim(tau$x_basis), c(n, length(tau$coefs)))
  expect_equal(dim(tau$x_basis_A0), c(n, length(tau$coefs)))
  expect_equal(dim(tau$x_basis_A1), c(n, length(tau$coefs)))
})

test_that("learn_tau works when external has controls only", {
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

  Pi <- learn_Pi(S, W, A, TRUE, "glm", 5, c(0.01, 0.99))
  theta <- learn_theta(W, A, Y, TRUE, "glm", 5, "gaussian", NULL)
  tau <- learn_tau(S, W, A, Y, Pi, theta, TRUE, "glmnet", 5, 1)

  expect_equal(length(tau$A0), n)
  expect_true(all(tau$A1 == 0))
  expect_equal(tau$pred, tau$A0) # ???
  expect_equal(dim(tau$x_basis), c(n, length(tau$coefs)))
  expect_equal(tau$x_basis, tau$x_basis_A0)
  expect_null(tau$x_basis_A1)
})
