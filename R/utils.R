learn_Q_S1 <- function(S, W, A, Y, method) {
  pred <- numeric(length = length(A))
  S1A1 <- numeric(length = length(A))
  S1A0 <- numeric(length = length(A))
  X <- data.frame(W, A = A)[S == 1,]
  X_S1A1 <- data.frame(W, A = 1)
  X_S1A0 <- data.frame(W, A = 0)

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = Y_bounded, keep = TRUE, alpha = 1, nfolds = 5, family = "gaussian")
    pred <- as.numeric(predict(fit, newx = as.matrix(X), s = "lambda.min", type = "response"))
    S1A1 <- as.numeric(predict(fit, newx = as.matrix(X_S1A1), s = "lambda.min", type = "response"))
    S1A0 <- as.numeric(predict(fit, newx = as.matrix(X_S1A0), s = "lambda.min", type = "response"))
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(X), Y = Y_bounded, family = "gaussian", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(X)))
    S1A1 <- as.numeric(predict(fit, new_data = as.matrix(X_S1A1)))
    S1A0 <- as.numeric(predict(fit, new_data = as.matrix(X_S1A0)))
  } else if (method == "glm") {
    fit <- glm(Y[S == 1] ~ ., data = X, family = "gaussian")
    pred <- as.numeric(predict(fit, newdata = data.frame(W, A = A), type = "response"))
    S1A1 <- as.numeric(predict(fit, newdata = X_S1A1, type = "response"))
    S1A0 <- as.numeric(predict(fit, newdata = X_S1A0, type = "response"))
  }

  return(list(pred = pred,
              S1A1 = S1A1,
              S1A0 = S1A0))
}

to_prob <- function(pred) {
  return(1 / (1 + exp(-pred)))
}

.bound <- function(X, bounds) {
  X[X < bounds[1]] <- bounds[1]
  X[X > bounds[2]] <- bounds[2]

  return(X)
}
