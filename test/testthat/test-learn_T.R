test_that("learn_T works", {
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

  g <- learn_g(S, W, A, 0.67, FALSE, "glm", 5, c(0.01, 0.99))
  theta_tilde <- learn_theta_tilde(W, Y, "glm", 5, "gaussian", NULL)
  T_est <- learn_T(W, A, Y, g, theta_tilde, "glmnet", 5)

  expect_equal(length(T_est$pred), n)
  expect_equal(dim(T_est$x_basis), c(n, length(T_est$coefs)))
})
