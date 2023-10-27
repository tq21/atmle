#' @title Learn nuisance function: trial enrollment mechanism
#'
#' @description Function to learn the trial enrollment mechanism
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
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
#' @importFrom purrr walk
#'
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param Pi_bounds A numeric vector of lower and upper bounds for the
#' trial enrollment probabilities. The first element is the lower bound,
#' and the second element is the upper bound.
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{The estimated trial enrollment probabilities;}
#' \item{A1}{The estimated trial enrollment probabilities under treatment;}
#' \item{A0}{The estimated trial enrollment probabilities under control.}
learn_Pi <- function(S,
                     W,
                     A,
                     controls_only,
                     method,
                     folds,
                     Pi_bounds,
                     cross_fit_nuisance = TRUE) {

  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- rep(NA, length(A))
  A1 <- rep(1, length(A))
  A0 <- rep(NA, length(A))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_Pi <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                            Lrnr_cv_selector$new(loss_loglik_binomial))
    task_Pi <- sl3_Task$new(data = data.table(W, A = A, S = S),
                            covariates = c(colnames(W), "A"),
                            outcome = "S", outcome_type = "binomial")
    fit_Pi <- lrnr_Pi$train(task_Pi)
    task_Pi_A0 <- sl3_Task$new(data = data.table(W, A = 0, S = S),
                               covariates = c(colnames(W), "A"),
                               outcome = "S", outcome_type = "binomial")
    A0 <- .bound(fit_Pi$predict(task_Pi_A0), Pi_bounds)
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1 # no treated in external
    } else {
      task_Pi_A1 <- sl3_Task$new(data = data.table(W, A = 1, S = S),
                                 covariates = c(colnames(W), "A"),
                                 outcome = "S", outcome_type = "binomial")
      A1 <- .bound(fit_Pi$predict(task_Pi_A1), Pi_bounds)
      pred[A == 1] <- A1[A == 1]
    }

  } else if (method == "glm") {

    if (controls_only) {
      # control
      X <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                      data = data.frame(W, A = A)))
      X_A0 <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                         data = data.frame(W, A = 0)))
      X_A1 <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                         data = data.frame(W, A = 1)))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- glm(S[train_idx] ~ ., data = X[train_idx,],
                     family = "binomial")
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit, newdata = X_A0[valid_idx,], type = "response")), Pi_bounds)
        })

        pred[A == 1] <- 1 # no treated in external
        A0 <- pred

      } else {
        # no cross fit
        fit <- glm(S ~ ., data = X, family = "binomial")
        A0 <- .bound(as.numeric(predict(fit, newdata = X_A0, type = "response")), Pi_bounds)
        pred <- A0
        pred[A == 1] <- 1 # no treated in external
      }

    } else {
      # treat + control
      X <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                      data = data.frame(W, A = A)))
      X_A0 <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                         data = data.frame(W, A = 0)))
      X_A1 <- as.data.frame(model.matrix(as.formula("~-1+.+A:."),
                                         data = data.frame(W, A = 1)))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- glm(S[train_idx] ~., data = X[train_idx, ],
                     family = "binomial")
          pred[valid_idx] <<- as.numeric(predict(
            fit, newdata = X[valid_idx,], type = "response"))
          A0[valid_idx] <<- as.numeric(predict(
            fit, newdata = X_A0[valid_idx,], type = "response"))
          A1[valid_idx] <<- as.numeric(predict(
            fit, newdata = X_A1[valid_idx,], type = "response"))
        })

      } else {
        # no cross fit
        fit <- glm(S ~ ., data = X, family = "binomial")
        A0 <- .bound(as.numeric(predict(fit, newdata = X_A0, type = "response")), Pi_bounds)
        A1 <- .bound(as.numeric(predict(fit, newdata = X_A1, type = "response")), Pi_bounds)
        pred <- .bound(as.numeric(predict(fit, newdata = X, type = "response")), Pi_bounds)
      }
    }

  } else if (method == "glmnet") {

    if (controls_only) {
      # control
      X <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A))
      X_A0 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 0))
      X_A1 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 1))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- cv.glmnet(x = as.matrix(X[train_idx, ]),
                           y = S[train_idx],
                           keep = TRUE, alpha = 1, nfolds = length(folds),
                           family = "binomial")
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit, newx = as.matrix(X_A0[valid_idx, ]), s = "lambda.min",
            type = "response")), Pi_bounds)
        })

        pred[A == 1] <- 1 # no treated in external
        A0 <- pred

      } else {
        # no cross fit
        fit <- cv.glmnet(x = X, y = S,
                         keep = TRUE, alpha = 1, nfolds = length(folds),
                         family = "binomial")
        A0 <- .bound(as.numeric(predict(
          fit, newx = X_A0, s = "lambda.min", type = "response")), Pi_bounds)
        pred <- A0
        pred[A == 1] <- 1 # no treated in external

      }

    } else {
      # treat + control
      X <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A))
      X_A0 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 0))
      X_A1 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 1))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- cv.glmnet(x = X[train_idx,], y = S[train_idx],
                           keep = TRUE, alpha = 1, nfolds = length(folds),
                           family = "binomial")
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit, newx = X[valid_idx,], s = "lambda.min", type = "response")), Pi_bounds)
          A0[valid_idx] <<- .bound(as.numeric(predict(
            fit, newx = X_A0[valid_idx,], s = "lambda.min", type = "response")), Pi_bounds)
          A1[valid_idx] <<- .bound(as.numeric(predict(
            fit, newx = X_A1[valid_idx,], s = "lambda.min", type = "response")), Pi_bounds)
        })

      } else {
        # no cross fit
        fit <- cv.glmnet(x = X, y = S, keep = TRUE, alpha = 1,
                         nfolds = length(folds), family = "binomial")
        A0 <- .bound(as.numeric(predict(
          fit, newx = X_A0, s = "lambda.min", type = "response")), Pi_bounds)
        A1 <- .bound(as.numeric(predict(
          fit, newx = X_A1, s = "lambda.min", type = "response")), Pi_bounds)
        pred <- .bound(as.numeric(predict(
          fit, newx = X, s = "lambda.min", type = "response")), Pi_bounds)
      }
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(list(pred = .bound(pred, Pi_bounds),
              A1 = .bound(A1, Pi_bounds),
              A0 = .bound(A0, Pi_bounds)))
}
