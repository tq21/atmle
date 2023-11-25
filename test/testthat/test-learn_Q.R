test_that("learn_Q works when outcome is continuous", {
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

  Q <- learn_Q(W, A, Y, "glm", 5, "gaussian", c(0, 0.5)) # glm
  expect_true(all(Q$pred >= 0 & Q$pred <= 0.5))
  expect_true(all(Q$A0 >= 0 & Q$A0 <= 0.5))
  expect_true(all(Q$A1 >= 0 & Q$A1 <= 0.5))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  Q <- learn_Q(W, A, Y, "glmnet", 5, "gaussian", c(-4, 4)) # glmnet
  expect_true(all(Q$pred >= -4 & Q$pred <= 4))
  expect_true(all(Q$A0 >= -4 & Q$A0 <= 4))
  expect_true(all(Q$A1 >= -4 & Q$A1 <= 4))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  Q <- learn_Q(W, A, Y, "sl3", 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(Q$pred >= -4 & Q$pred <= 4))
  expect_true(all(Q$A0 >= -4 & Q$A0 <= 4))
  expect_true(all(Q$A1 >= -4 & Q$A1 <= 4))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  Q <- learn_Q(W, A, Y, lrnrs, 5, "gaussian", c(-4, 4)) # sl3
  expect_true(all(Q$pred >= -4 & Q$pred <= 4))
  expect_true(all(Q$A0 >= -4 & Q$A0 <= 4))
  expect_true(all(Q$A1 >= -4 & Q$A1 <= 4))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  expect_error(learn_Q(W, A, Y, "abc", 5, "gaussian", c(-4, 4)))
})

test_that("learn_Q works when outcome is binary", {
  # simulate data
  set.seed(123)
  n <- 500
  S <- rbinom(n, 1, 0.5)
  W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.67)
  A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
  Y <- rbinom(n, 1, plogis(-0.5-0.8*W1-1.1*W2+1.5*A+(1-S)*(0.9+2.6*W1)))

  Q <- learn_Q(W, A, Y, "glm", 5, "binomial", c(0.2, 0.8)) # glm
  expect_true(all(Q$pred >= 0.2 & Q$pred <= 0.8))
  expect_true(all(Q$A0 >= 0.2 & Q$A0 <= 0.8))
  expect_true(all(Q$A1 >= 0.2 & Q$A1 <= 0.8))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  Q <- learn_Q(W, A, Y, "glmnet", 5, "binomial", NULL) # glmnet
  expect_true(all(Q$pred >= 0 & Q$pred <= 1))
  expect_true(all(Q$A0 >= 0 & Q$A0 <= 1))
  expect_true(all(Q$A1 >= 0 & Q$A1 <= 1))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  Q <- learn_Q(W, A, Y, "sl3", 5, "binomial", NULL) # sl3
  expect_true(all(Q$pred >= 0 & Q$pred <= 1))
  expect_true(all(Q$A0 >= 0 & Q$A0 <= 1))
  expect_true(all(Q$A1 >= 0 & Q$A1 <= 1))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  lrnrs <- list(Lrnr_mean$new(), Lrnr_glm$new())
  Q <- learn_Q(W, A, Y, lrnrs, 5, "binomial", NULL) # sl3
  expect_true(all(Q$pred >= 0 & Q$pred <= 1))
  expect_true(all(Q$A0 >= 0 & Q$A0 <= 1))
  expect_true(all(Q$A1 >= 0 & Q$A1 <= 1))
  expect_equal(Q$pred[A == 0], Q$A0[A == 0])
  expect_equal(Q$pred[A == 1], Q$A1[A == 1])

  expect_error(learn_Q(W, A, Y, "abc", 5, "binomial", NULL))
})
