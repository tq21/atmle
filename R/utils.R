`%+%` <- function(a, b) paste0(a, b)

svd_pseudo_inv <- function(mat, tol = 1e-3) {
  svd_res <- svd(mat)
  D_inv <- ifelse(svd_res$d > tol, 1 / svd_res$d, 0)
  return(svd_res$v %*% diag(D_inv) %*% t(svd_res$u))
}

learn_procova_score <- function(S,
                                W,
                                A,
                                Y,
                                delta,
                                family,
                                controls_only,
                                v_folds,
                                method) {
  if (method == "glmnet") {
    if (controls_only) {
      X_train <- as.matrix(W[S == 0 & delta == 1, , drop = FALSE])
      X_eval <- as.matrix(W[S == 1, , drop = FALSE])
      fit <- cv.glmnet(
        x = X_train,
        y = Y[S == 0 & delta == 1],
        keep = TRUE, alpha = 1, nfolds = v_folds, family = family
      )
      pred <- as.numeric(predict(
        fit,
        newx = X_eval, s = "lambda.min", type = "response"
      ))
    } else {
      X_train <- as.matrix(data.frame(W, A = A)[S == 0 & delta == 1, ])
      X_eval <- as.matrix(data.frame(W, A = A)[S == 1, ])
      fit <- cv.glmnet(
        x = X_train,
        y = Y[S == 0 & delta == 1],
        keep = TRUE, alpha = 1, nfolds = v_folds, family = family
      )
      pred <- as.numeric(predict(
        fit,
        newx = X_eval, s = "lambda.min", type = "response"
      ))
    }
  } else if (method == "glm") {
    if (controls_only) {
      X_train <- W[S == 0 & delta == 1, , drop = FALSE]
      X_eval <- W[S == 1, , drop = FALSE]
      fit <- glm(Y[S == 0 & delta == 1] ~ ., data = X_train, family = family)
      pred <- as.numeric(predict(fit, newdata = X_eval, type = "response"))
    } else {
      X_train <- data.frame(W, A = A)[S == 0 & delta == 1, ]
      X_eval <- data.frame(W, A = A)[S == 1, ]
      fit <- glm(Y[S == 0 & delta == 1] ~ ., data = X_train, family = family)
      pred <- as.numeric(predict(fit, newdata = X_eval, type = "response"))
    }
  }

  return(pred)
}

to_prob <- function(pred) {
  return(1 / (1 + exp(-pred)))
}

bound <- function(X) {
  X_max <- max(X, na.rm = TRUE)
  X_min <- min(X, na.rm = TRUE)

  return((X - X_min) / (X_max - X_min))
}
