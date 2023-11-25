test_that("basic functionality test, external has both treated and controls", {
  # simulate data
  set.seed(123)
  n <- 2000
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
  data <- data.frame(S, W1, W2, A, Y)

  # test atmle_pooled = TRUE
  res <- atmle(data,
               S_node = c(1),
               W_node = c(2, 3),
               A_node = 4,
               Y_node = 5,
               controls_only = FALSE,
               family = "gaussian",
               atmle_pooled = TRUE,
               g_rct = 0.67,
               verbose = FALSE)
  expect_equal(res$est, 1.5, tolerance = 1)

  # test atmle_pooled = FALSE
  res <- atmle(data,
               S_node = c(1),
               W_node = c(2, 3),
               A_node = 4,
               Y_node = 5,
               controls_only = FALSE,
               family = "gaussian",
               atmle_pooled = FALSE,
               g_rct = 0.67,
               verbose = FALSE)
  expect_equal(res$est, 1.5, tolerance = 1)
})

test_that("basic functionality test, external has only controls", {
  # simulate data
  set.seed(123)
  n <- 2000
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
  data <- data.frame(S, W1, W2, A, Y)

  # test atmle_pooled = TRUE
  res <- atmle(data,
               S_node = c(1),
               W_node = c(2, 3),
               A_node = 4,
               Y_node = 5,
               controls_only = TRUE,
               family = "gaussian",
               atmle_pooled = TRUE,
               g_rct = 0.67,
               verbose = FALSE)
  expect_equal(res$est, 1.5, tolerance = 1)

  # test atmle_pooled = FALSE
  res <- atmle(data,
               S_node = c(1),
               W_node = c(2, 3),
               A_node = 4,
               Y_node = 5,
               controls_only = TRUE,
               family = "gaussian",
               atmle_pooled = FALSE,
               g_rct = 0.67,
               verbose = FALSE)
  expect_equal(res$est, 1.5, tolerance = 1)
})
