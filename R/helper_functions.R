library(hal9001)
library(data.table)
library(tmle)
library(glmnet)
library(origami)
library(purrr)

Q_tmle <- function(g, Q, A, Y_bound) {
  wt <- A / g + (1 - A) / (1 - g)
  H1W <- A
  H0W <- 1 - A

  submodel <- glm(Y_bound ~ -1 + offset(Q$pred) + H0W + H1W, family = "quasibinomial", weights = wt)
  epsilon <- coef(submodel)

  Q_star <- Q$pred + epsilon[1] * H0W + epsilon[2] * H1W
  Q_A1_star <- Q$A1 + rep(epsilon[1], length(Y_bound))
  Q_A0_star <- Q$A0 + rep(epsilon[2], length(Y_bound))

  return(list(
    Q_star = Q_star,
    A1 = Q_A1_star,
    A0 = Q_A0_star
  ))
}

.bound <- function(x, bounds) {
  return(pmin(pmax(x, bounds[1]), bounds[2]))
}

learn_g_S1 <- function(S, W, A, g_rct, method = "glmnet") {
  pred <- numeric(length = length(A))
  pred[S == 1] <- p_rct
  X <- data.frame(W[S == 0, ])

  if (method == "lasso") {
    fit <- cv.glmnet(x = as.matrix(X), y = A[S == 0], keep = TRUE, alpha = 1, nfolds = 5, family = "binomial")
    pred[S == 0] <- as.numeric(predict(fit, newx = as.matrix(X), s = "lambda.min", type = "response"))
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(X), Y = A[S == 0], family = "binomial", smoothness_orders = 0)
    pred[S == 0] <- as.numeric(predict(fit, new_data = as.matrix(X), type = "response"))
  } else if (method == "glm") {
    fit <- glm(A[S == 0] ~ ., data = X, family = "binomial")
    pred[S == 0] <- as.numeric(predict(fit, newdata = X, type = "response"))
  }

  return(list(pred = pred))
}

learn_S_W <- function(S, W, method = "glmnet") {
  pred <- numeric(length = length(S))
  X <- data.frame(W)

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(X), y = S, keep = TRUE, alpha = 1, nfolds = 5, family = "binomial")
    pred <- as.numeric(predict(fit, newx = as.matrix(X), s = "lambda.min", type = "response"))
  } else if (method == "HAL") {
    fit <- fit_hal(X = as.matrix(X), Y = S, family = "binomial", smoothness_orders = 0)
    pred <- as.numeric(predict(fit, new_data = as.matrix(X), type = "response"))
  } else if (method == "glm") {
    fit <- glm(S ~ ., data = X, family = "binomial")
    pred <- as.numeric(predict(fit, newdata = X, type = "response"))
  }

  return(list(pred = pred))
}

target_Q <- function(S, W, A, Y, Pi, g, Q, delta, g_delta) {
  # bound Y
  min_Y <- min(Y, Q$pred, Q$S1A1, Q$S1A0, na.rm = TRUE) - 0.001
  max_Y <- max(Y, Q$pred, Q$S1A1, Q$S1A0, na.rm = TRUE) + 0.001
  Y_bounded <- (Y - min_Y) / (max_Y - min_Y)
  Q$pred <- (Q$pred - min_Y) / (max_Y - min_Y)
  Q$S1A1 <- (Q$S1A1 - min_Y) / (max_Y - min_Y)
  Q$S1A0 <- (Q$S1A0 - min_Y) / (max_Y - min_Y)

  # clever covariates
  wt <- S / Pi$pred
  H1_n <- wt * (A / g) * (delta / g_delta$pred)
  H0_n <- wt * ((1 - A) / (1 - g)) * (delta / g_delta$pred)

  # logistic submodel
  epsilon <- coef(glm(Y_bounded ~ -1 + offset(qlogis(Q$pred)) + H0_n + H1_n, family = "quasibinomial"))
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Q_star <- NULL
  Q_star$pred <- plogis(qlogis(Q$pred) + epsilon[1] * H0_n + epsilon[2] * H1_n)
  Q_star$S1A0 <- plogis(qlogis(Q$S1A0) + epsilon[1] * H0_n)
  Q_star$S1A1 <- plogis(qlogis(Q$S1A1) + epsilon[2] * H1_n)

  # scale back
  Q_star$pred <- Q_star$pred * (max_Y - min_Y) + min_Y
  Q_star$S1A0 <- Q_star$S1A0 * (max_Y - min_Y) + min_Y
  Q_star$S1A1 <- Q_star$S1A1 * (max_Y - min_Y) + min_Y

  return(Q_star)
}
