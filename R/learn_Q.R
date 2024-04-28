#' @title Learn outcome regression
#'
#' @description Function to learn the outcome regression
#' \eqn{Q(A,W)=\mathbb{E}(Y\mid A,W)}.
#'
#' @keywords nuisance
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table
#' @importFrom sl3 Stack
#' @importFrom sl3 make_learner
#' @importFrom sl3 sl3_Task
#' @importFrom sl3 Pipeline
#' @importFrom sl3 Lrnr_cv
#' @importFrom sl3 Lrnr_cv_selector
#' @importFrom sl3 loss_loglik_binomial
#'
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Either \code{"gaussian"} or \code{"binomial"}.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates and treatment.
#' The first element is the lower bound, and the second element is the upper
#' bound.
#' @returns A numeric vector of estimated values.
learn_Q <- function(W,
                    A,
                    Y,
                    delta,
                    method,
                    v_folds,
                    family,
                    theta_bounds) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  if (is.null(theta_bounds)) {
    if (family == "gaussian") {
      theta_bounds <- c(-Inf, Inf)
    } else if (family == "binomial") {
      theta_bounds <- c(0.01, 0.99)
    }
  }

  pred <- numeric(length(A))
  A0 <- numeric(length = length(A))
  A1 <- numeric(length = length(A))
  X <- data.frame(W, A = A)
  X_A0 <- data.frame(W, A = 0)
  X_A1 <- data.frame(W, A = 1)

  if (is.list(method)) {
    if (family == "gaussian") {
      lrnr_stack <- Stack$new(method)
      lrnr_Q <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_squared_error)
      )
      task_Q <- sl3_Task$new(
        data = data.table(W, Y = Y, A = A)[delta == 1],
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "continuous"
      )
      fit_Q <- lrnr_Q$train(task_Q)

      task_Q_A0 <- sl3_Task$new(
        data = data.table(W, Y = Y, A = 0),
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "continuous"
      )
      A0 <- .bound(fit_Q$predict(task_Q_A0), theta_bounds)
      pred[A == 0] <- A0[A == 0]
      task_Q_A1 <- sl3_Task$new(
        data = data.table(W, Y = Y, A = 1),
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "continuous"
      )
      A1 <- .bound(fit_Q$predict(task_Q_A1), theta_bounds)
      pred[A == 1] <- A1[A == 1]
    } else if (family == "binomial") {
      lrnr_stack <- Stack$new(method)
      lrnr_Q <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_loglik_binomial)
      )
      task_Q <- sl3_Task$new(
        data = data.table(W, Y = Y, A = A)[delta == 1],
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "binomial"
      )
      fit_Q <- lrnr_Q$train(task_Q)

      task_Q_A0 <- sl3_Task$new(
        data = data.table(W, Y = Y, A = 0),
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "binomial"
      )
      A0 <- .bound(fit_Q$predict(task_Q_A0), theta_bounds)
      pred[A == 0] <- A0[A == 0]
      task_Q_A1 <- sl3_Task$new(
        data = data.table(W, Y = Y, A = 1),
        covariates = c(colnames(W), "A"),
        outcome = "Y", outcome_type = "binomial"
      )
      A1 <- .bound(fit_Q$predict(task_Q_A1), theta_bounds)
      pred[A == 1] <- A1[A == 1]
    } else {
      stop("Invalid family. Must be either 'gaussian' or 'binomial'.")
    }
  } else if (method == "glm") {
    fit <- glm(Y[delta == 1] ~ ., data = X[delta == 1, , drop = FALSE], family = family)
    A0 <- .bound(
      as.numeric(predict(fit, newdata = X_A0, type = "response")),
      theta_bounds
    )
    pred[A == 0] <- A0[A == 0]
    A1 <- .bound(
      as.numeric(predict(fit, newdata = X_A1, type = "response")),
      theta_bounds
    )
    pred[A == 1] <- A1[A == 1]
  } else if (method == "glmnet") {
    fit <- cv.glmnet(
      x = as.matrix(X[delta == 1, , drop = FALSE]), y = Y[delta == 1],
      keep = TRUE, alpha = 1, nfolds = v_folds,
      family = family
    )
    A0 <- .bound(
      as.numeric(predict(fit,
        newx = as.matrix(X_A0),
        s = "lambda.min", type = "response"
      )),
      theta_bounds
    )
    pred[A == 0] <- A0[A == 0]
    A1 <- .bound(
      as.numeric(predict(fit,
        newx = as.matrix(X_A1),
        s = "lambda.min", type = "response"
      )),
      theta_bounds
    )
    pred[A == 1] <- A1[A == 1]
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(list(
    pred = pred,
    A1 = A1,
    A0 = A0
  ))
}
