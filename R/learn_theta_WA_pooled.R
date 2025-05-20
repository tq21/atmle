#' @title Learn nuisance function: conditional mean of outcome given baseline
#' covariates
#'
#' @description Function to learn the conditional mean of outcome given
#' baseline covariates, \eqn{\tilde{\theta}(W)=\mathbb{E}(Y\mid W)}.
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
#' @importFrom sl3 loss_squared_error
#' @importFrom purrr walk
#'
#' @param W A matrix of baseline covariates.
#' @param Y A vector of outcomes.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Either \code{"gaussian"} or \code{"binomial"}.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates.
#' The first element is the lower bound, and the second element is the upper
#' bound.
#'
#' @returns A numeric vector of the estimated values.
learn_theta_WA_pooled <- function(W,
                                  Y,
                                  delta,
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

  pred <- numeric(length(Y))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_theta_tilde <- NULL
    task_theta_tilde <- NULL
    if (family == "gaussian") {
      lrnr_theta_tilde <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_squared_error)
      )
      task_train <- sl3_Task$new(
        data = data.table(W, Y = Y)[delta == 1],
        covariates = colnames(W), outcome = "Y", outcome_type = "continuous"
      )
      Y_tmp <- Y
      Y_tmp[delta == 0] <- 0
      task_pred <- sl3_Task$new(
        data = data.table(W, Y = Y_tmp),
        covariates = colnames(W), outcome = "Y", outcome_type = "continuous"
      )
    } else if (family == "binomial") {
      lrnr_theta_tilde <- make_learner(
        Pipeline, Lrnr_cv$new(lrnr_stack),
        Lrnr_cv_selector$new(loss_loglik_binomial)
      )
      task_train <- sl3_Task$new(
        data = data.table(W, Y = Y)[delta == 1],
        covariates = colnames(W), outcome = "Y", outcome_type = "binomial"
      )
      Y_tmp <- Y
      Y_tmp[delta == 0] <- 0
      task_pred <- sl3_Task$new(
        data = data.table(W, Y = Y_tmp),
        covariates = colnames(W), outcome = "Y", outcome_type = "binomial"
      )
    }

    fit_theta_tilde <- lrnr_theta_tilde$train(task_train)
    pred <- .bound(fit_theta_tilde$predict(task_pred), theta_bounds)
  } else if (method == "glm") {
    X <- data.frame(W) # USE MODEL MATRIX, SOMETIMES CHARACTERS MESS UP

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
      pred <- .bound(
        as.numeric(predict(fit, newdata = X, type = "response")),
        theta_bounds
      )
    }
  } else if (method == "glmnet") {
    X <- as.matrix(W)

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(
          x = X[train_idx, ][delta[train_idx] == 1, ],
          y = Y[train_idx][delta[train_idx] == 1], keep = TRUE,
          alpha = 1, nfolds = length(folds), family = family
        )
        pred[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newx = X[valid_idx, ], s = "lambda.min", type = "response"
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
        newx = X, s = "lambda.min", type = "response"
      )), theta_bounds)
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(pred)
}
