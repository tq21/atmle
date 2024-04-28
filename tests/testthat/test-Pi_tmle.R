test_that("Pi_tmle works when external has both treated and controls", {
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
  Y <- -0.5 - 0.8 * W1 - 1.1 * W2 + 1.5 * A + UY + (1 - S) * (0.9 + 2.6 * W1)
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # estimate nuisance parameters
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = FALSE,
    method = "glm",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  )
  Pi <- learn_Pi(
    S = S,
    W = W,
    A = A,
    controls_only = FALSE,
    method = "glm",
    folds = folds,
    Pi_bounds = c(0.01, 0.99)
  )
  theta <- learn_theta(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    controls_only = FALSE,
    method = "glm",
    folds = folds,
    family = "gaussian",
    theta_bounds = NULL
  )
  tau <- learn_tau(
    S = S,
    W = W,
    A = A,
    Y = Y,
    Pi = Pi,
    theta = theta,
    g = g,
    delta = delta,
    controls_only = FALSE,
    method = "glmnet",
    v_folds = 5,
    max_degree = 1,
    min_working_model = FALSE,
    min_working_model_screen = FALSE,
    undersmooth = 0,
    target_gwt = TRUE,
    Pi_bounds = c(0.01, 0.99),
    weights = delta
  )

  # target_gwt = TRUE
  Pi_star <- Pi_tmle(
    S = S,
    W = W,
    A = A,
    g = g,
    tau = tau,
    Pi = Pi,
    controls_only = FALSE,
    target_gwt = TRUE,
    Pi_bounds = c(0.2, 0.8)
  )
  expect_length(Pi_star$pred, n)
  expect_length(Pi_star$A0, n)
  expect_length(Pi_star$A1, n)
  expect_true(all(Pi_star$pred >= 0.2 & Pi_star$pred <= 0.8))
  expect_true(all(Pi_star$A0 >= 0.2 & Pi_star$A0 <= 0.8))
  expect_true(all(Pi_star$A1 >= 0.2 & Pi_star$A1 <= 0.8))

  # target_gwt = FALSE
  Pi_star <- Pi_tmle(
    S = S,
    W = W,
    A = A,
    g = g,
    tau = tau,
    Pi = Pi,
    controls_only = FALSE,
    target_gwt = FALSE,
    Pi_bounds = c(0.2, 0.8)
  )
  expect_length(Pi_star$pred, n)
  expect_length(Pi_star$A0, n)
  expect_length(Pi_star$A1, n)
  expect_true(all(Pi_star$pred >= 0.2 & Pi_star$pred <= 0.8))
  expect_true(all(Pi_star$A0 >= 0.2 & Pi_star$A0 <= 0.8))
  expect_true(all(Pi_star$A1 >= 0.2 & Pi_star$A1 <= 0.8))
})

test_that("Pi_tmle works when external has controls only", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0
  UY <- rnorm(n, 0, 1)
  U_bias <- rnorm(n, 0, 0.5)
  Y <- -0.5 - 0.8 * W1 - 1.1 * W2 + 1.5 * A + UY + (1 - S) * (0.9 + 2.6 * W1)
  delta <- rep(1, n)

  # make cv folds
  cv_strata <- paste0(S, "-", A)
  suppressWarnings({
    folds <- make_folds(
      n = n, V = 5,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # estimate nuisance parameters
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    g_rct = 0.67,
    controls_only = TRUE,
    method = "glm",
    folds = folds,
    g_bounds = c(0.01, 0.99)
  )
  Pi <- learn_Pi(
    S = S,
    W = W,
    A = A,
    controls_only = TRUE,
    method = "glm",
    folds = folds,
    Pi_bounds = c(0.01, 0.99)
  )
  theta <- learn_theta(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    controls_only = TRUE,
    method = "glm",
    folds = folds,
    family = "gaussian",
    theta_bounds = NULL
  )
  tau <- learn_tau(
    S = S,
    W = W,
    A = A,
    Y = Y,
    Pi = Pi,
    theta = theta,
    g = g,
    delta = delta,
    controls_only = TRUE,
    method = "glmnet",
    v_folds = 5,
    max_degree = 1,
    min_working_model = FALSE,
    min_working_model_screen = FALSE,
    undersmooth = 0,
    target_gwt = TRUE,
    Pi_bounds = c(0.01, 0.99),
    weights = delta
  )

  # target_gwt = TRUE
  Pi_star <- Pi_tmle(
    S = S,
    W = W,
    A = A,
    g = g,
    tau = tau,
    Pi = Pi,
    controls_only = TRUE,
    target_gwt = TRUE,
    Pi_bounds = c(0.2, 0.8)
  )
  expect_length(Pi_star$pred, n)
  expect_length(Pi_star$A0, n)
  expect_length(Pi_star$A1, n)
  expect_true(all(Pi_star$A1 == 1))
  expect_true(all(Pi_star$pred[A == 1] == 1))
  expect_true(all(Pi_star$pred[A == 0] >= 0.2 & Pi_star$pred[A == 0] <= 0.8))

  Pi_star <- Pi_tmle(
    S = S,
    W = W,
    A = A,
    g = g,
    tau = tau,
    Pi = Pi,
    controls_only = TRUE,
    target_gwt = FALSE,
    Pi_bounds = c(0.2, 0.8)
  )
  expect_length(Pi_star$pred, n)
  expect_length(Pi_star$A0, n)
  expect_length(Pi_star$A1, n)
  expect_true(all(Pi_star$A1 == 1))
  expect_true(all(Pi_star$pred[A == 1] == 1))
  expect_true(all(Pi_star$pred[A == 0] >= 0.2 & Pi_star$pred[A == 0] <= 0.8))
})
