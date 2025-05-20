#' @title Learn nuisance function: conditional mean of outcome given baseline
#' covariates and treatment
#'
#' @description Function to learn the conditional mean of outcome given
#' baseline covariates and treatment,
#' \eqn{\theta(W,A)=\mathbb{E}(Y\mid W,A)}.
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
#' @importFrom sl3 loss_squared_error
#' @importFrom purrr walk
#'
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators.
#' @param Y A vector of outcomes.
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
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
#'
#' @export
#'
#' @returns A numeric vector of the estimated values.
learn_theta_WA <- function(W,
                           A,
                           Y,
                           delta,
                           controls_only,
                           method,
                           folds,
                           family,
                           theta_bounds,
                           cross_fit_nuisance) {
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

  pred <- rep(NA, length(A))

  if (is.list(method)) {
    if (family == "gaussian") {
      lrnr_stack <- Stack$new(method)
      lrnr_theta <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_squared_error)
      )

      if (controls_only) {
        task_train <- sl3_Task$new(
          data = data.table(W, Y = Y)[A == 0 & delta == 1],
          covariates = colnames(W), outcome = "Y", outcome_type = "continuous"
        )
        Y_tmp <- Y
        Y_tmp[delta == 0] <- 0
        task_pred <- sl3_Task$new(
          data = data.table(W, Y = Y_tmp)[A == 0],
          covariates = colnames(W), outcome = "Y", outcome_type = "continuous"
        )
        fit_theta <- lrnr_theta$train(task_train)
        pred[A == 0] <- .bound(fit_theta$predict(task_pred), theta_bounds)
      } else {
        task_train <- sl3_Task$new(
          data = data.table(W, Y = Y, A = A)[delta == 1],
          covariates = c(colnames(W), "A"),
          outcome = "Y", outcome_type = "continuous"
        )
        Y_tmp <- Y
        Y_tmp[delta == 0] <- 0
        task_pred <- sl3_Task$new(
          data = data.table(W, Y = Y_tmp, A = A),
          covariates = c(colnames(W), "A"),
          outcome = "Y", outcome_type = "continuous"
        )
        fit_theta <- lrnr_theta$train(task_train)
        pred <- .bound(fit_theta$predict(task_pred), theta_bounds)
      }
    } else if (family == "binomial") {
      lrnr_stack <- Stack$new(method)
      lrnr_theta <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_loglik_binomial)
      )

      if (controls_only) {
        task_train <- sl3_Task$new(
          data = data.table(W, Y = Y)[A == 0 & delta == 1],
          covariates = colnames(W), outcome = "Y", outcome_type = "binomial"
        )
        Y_tmp <- Y
        Y_tmp[delta == 0] <- 0
        task_pred <- sl3_Task$new(
          data = data.table(W, Y = Y_tmp)[A == 0],
          covariates = colnames(W), outcome = "Y", outcome_type = "binomial"
        )
        fit_theta <- lrnr_theta$train(task_train)
        pred[A == 0] <- .bound(fit_theta$predict(task_pred), theta_bounds)
      } else {
        task_train <- sl3_Task$new(
          data = data.table(W, Y = Y, A = A)[delta == 1],
          covariates = c(colnames(W), "A"),
          outcome = "Y", outcome_type = "binomial"
        )
        Y_tmp <- Y
        Y_tmp[delta == 0] <- 0
        task_pred <- sl3_Task$new(
          data = data.table(W, Y = Y_tmp, A = A),
          covariates = c(colnames(W), "A"),
          outcome = "Y", outcome_type = "binomial"
        )
        fit_theta <- lrnr_theta$train(task_train)
        pred <- .bound(fit_theta$predict(task_pred), theta_bounds)
      }
    } else {
      stop("Invalid family. Must be either 'gaussian' or 'binomial'.")
    }
  } else if (method == "glm") {
    if (controls_only) {
      # control
      X <- as.data.frame(model.matrix(
        as.formula("~-1+.+A:."),
        data = data.frame(W, A = A)
      ))
      X_A0 <- as.data.frame(model.matrix(
        as.formula("~-1+.+A:."),
        data = data.frame(W, A = 0)
      ))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- glm(Y[train_idx][delta[train_idx] == 1] ~ .,
            data = X[train_idx, ][delta[train_idx] == 1, ],
            family = family
          )
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit,
            newdata = X_A0[valid_idx, ], type = "response"
          )), theta_bounds)
        })
      } else {
        # no cross fit
        fit <- glm(Y[delta == 1] ~ ., data = X[delta == 1, ], family = family)
        pred <- .bound(as.numeric(predict(fit, newdata = X_A0, type = "response")), theta_bounds)
      }
    } else {
      # treat + control
      X <- as.data.frame(model.matrix(
        as.formula("~-1+.+A:."),
        data = data.frame(W, A = A)
      ))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- glm(Y[train_idx][delta[train_idx] == 1] ~ .,
            data = X[train_idx, ][delta[train_idx] == 1, ],
            family = family
          )
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit,
            newdata = X[valid_idx, ], type = "response"
          )), theta_bounds)
        })
      } else {
        # no cross fit
        fit <- glm(Y[delta == 1] ~ ., data = X[delta == 1, ], family = family)
        pred <- .bound(as.numeric(predict(fit, newdata = X, type = "response")), theta_bounds)
      }
    }
  } else if (method == "glmnet") {
    if (controls_only) {
      # control
      X <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A))
      X_A0 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 0))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- cv.glmnet(
            x = X[train_idx, ][delta[train_idx] == 1, ],
            y = Y[train_idx][delta[train_idx] == 1],
            keep = TRUE, alpha = 1, nfolds = length(folds),
            family = family
          )
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit,
            newx = X_A0[valid_idx, ], s = "lambda.min", type = "response"
          )), theta_bounds)
        })
      } else {
        # no cross fit
        fit <- cv.glmnet(
          x = X[delta == 1, ], y = Y[delta == 1], keep = TRUE,
          alpha = 1, nfolds = length(folds), family = family
        )
        pred <- .bound(as.numeric(predict(
          fit,
          newx = X_A0, s = "lambda.min", type = "response"
        )), theta_bounds)
      }
    } else {
      # treat + control
      X <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit <- cv.glmnet(
            x = X[train_idx, ][delta[train_idx] == 1, ],
            y = Y[train_idx][delta[train_idx] == 1],
            keep = TRUE, alpha = 1, nfolds = length(folds),
            family = family
          )
          pred[valid_idx] <<- .bound(as.numeric(predict(
            fit,
            newx = X[valid_idx, ], s = "lambda.min", type = "response"
          )), theta_bounds)
        })
      } else {
        # no cross fit
        fit <- cv.glmnet(
          x = X[delta == 1, ], y = Y[delta == 1],
          keep = TRUE, alpha = 1, nfolds = length(folds),
          family = family
        )
        pred <- .bound(as.numeric(predict(
          fit,
          newx = X, s = "lambda.min", type = "response"
        )), theta_bounds)
      }
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(pred)
}
