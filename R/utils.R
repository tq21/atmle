
`%+%` <- function(a, b) paste0(a, b)

learn_Q_S1 <- function(S, W, A, Y, delta, method = "glm") {
  pred <- numeric(length = length(A))
  S1A1 <- numeric(length = length(A))
  S1A0 <- numeric(length = length(A))
  X <- data.frame(W, A = A)[S == 1 & delta == 1,]
  X_S1A1 <- data.frame(W, A = 1)
  X_S1A0 <- data.frame(W, A = 0)

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = Y[S == 1 & delta == 1],
                     keep = TRUE, alpha = 1, nfolds = 5, family = "gaussian")
    pred <- as.numeric(predict(fit, newx = as.matrix(X), s = "lambda.min", type = "response"))
    S1A1 <- as.numeric(predict(fit, newx = as.matrix(X_S1A1), s = "lambda.min", type = "response"))
    S1A0 <- as.numeric(predict(fit, newx = as.matrix(X_S1A0), s = "lambda.min", type = "response"))
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(X), Y = Y[delta == 1],
                   family = "gaussian", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(X)))
    S1A1 <- as.numeric(predict(fit, new_data = as.matrix(X_S1A1)))
    S1A0 <- as.numeric(predict(fit, new_data = as.matrix(X_S1A0)))
  } else if (method == "glm") {
    fit <- glm(Y[S == 1 & delta == 1] ~ ., data = X, family = "gaussian")
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

bound <- function(X) {
  X_max <- max(X, na.rm = TRUE)
  X_min <- min(X, na.rm = TRUE)

  return((X-X_min)/(X_max-X_min))
}
