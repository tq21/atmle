test_that("learn_g works when external has both treated and controls", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))

  g <- learn_g(S, W, A, 0.67, FALSE, "glm", 5, c(0.01, 0.5)) # glm
  expect_true(min(g) >= 0.01 & max(g) <= 0.5)

  g <- learn_g(S, W, A, 0.67, FALSE, "glmnet", 5, c(0.01, 0.99)) # glmnet
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)

  g <- learn_g(S, W, A, 0.67, FALSE, "sl3", 5, c(0.01, 0.99)) # sl3
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  g <- learn_g(S, W, A, 0.67, FALSE, lrnrs, 5, c(0.01, 0.99)) # sl3 learners
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)

  expect_error(learn_g(S, W, A, 0.67, FALSE, "abc", 5, c(0.01, 0.99)))
})

test_that("learn_g works when external has controls only", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0

  g <- learn_g(S, W, A, 0.67, TRUE, "glm", 5, c(0.01, 0.05)) # glm
  expect_true(min(g) >= 0.01 & max(g) <= 0.05)

  g <- learn_g(S, W, A, 0.67, TRUE, "glmnet", 5, c(0.01, 0.05)) # glmnet
  expect_true(min(g) >= 0.01 & max(g) <= 0.05)

  g <- learn_g(S, W, A, 0.67, TRUE, "sl3", 5, c(0.01, 0.99)) # sl3
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  g <- learn_g(S, W, A, 0.67, TRUE, lrnrs, 5, c(0.01, 0.99)) # sl3 learners
  expect_true(min(g) >= 0.01 & max(g) <= 0.99)
})
