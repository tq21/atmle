test_that("learn_Pi works when external has both treated and controls", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))

  Pi <- learn_Pi(S, W, A, FALSE, "glm", 5, c(0.01, 0.2)) # glm
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.2))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.2))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.2))

  Pi <- learn_Pi(S, W, A, FALSE, "glmnet", 5, c(0.01, 0.99)) # glmnet
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))

  Pi <- learn_Pi(S, W, A, FALSE, "sl3", 5, c(0.01, 0.99)) # sl3
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  Pi <- learn_Pi(S, W, A, FALSE, lrnrs, 5, c(0.01, 0.99)) # sl3 learners
  expect_true(all(Pi$pred >= 0.01 & Pi$pred <= 0.99))
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$A1 >= 0.01 & Pi$A1 <= 0.99))

  expect_error(learn_Pi(S, W, A, FALSE, "abc", 5, c(0.01, 0.99)))
})

test_that("learn_Pi works when external has controls only", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- 0

  Pi <- learn_Pi(S, W, A, TRUE, "glm", 5, c(0.01, 0.2)) # glm
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.2))
  expect_true(all(Pi$pred[A == 1] == 1))

  Pi <- learn_Pi(S, W, A, TRUE, "glmnet", 5, c(0.01, 0.99)) # glmnet
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$pred[A == 1] == 1))

  Pi <- learn_Pi(S, W, A, TRUE, "sl3", 5, c(0.01, 0.99)) # sl3
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$pred[A == 1] == 1))

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  Pi <- learn_Pi(S, W, A, TRUE, lrnrs, 5, c(0.01, 0.99)) # sl3 learners
  expect_true(all(Pi$A0 >= 0.01 & Pi$A0 <= 0.99))
  expect_true(all(Pi$pred[A == 1] == 1))

  expect_error(learn_Pi(S, W, A, TRUE, "abc", 5, c(0.01, 0.99)))
})
