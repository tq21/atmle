#' @title Learn nuisance function: treatment mechanism
#'
#' @description Function to learn the treatment mechanism
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
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
#' @param g_rct A numeric of treatment probability in RCT.
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param g_bounds A numeric vector of lower and upper bounds for the
#' treatment mechanism. The first element is the lower bound, and the second
#' element is the upper bound.
#'
#' @returns A numeric vector of estimated treatment probabilities.
#'
#' @details We provide some details on estimating the treatment mechanism
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}. Note that, we can rewrite the
#' treatment mechanism into the following form:
#' \deqn{g(A\mid W)=\mathbb{P}(A=1\mid S=1,W)\mathbb{P}(S=1\mid W)+
#' \mathbb{P}(A=1\mid S=0,W)\mathbb{P}(S=0\mid W).}
#' This form allows us to use the RCT randomization probability
#' (typically known) for the term \eqn{\mathbb{P}(A=1\mid S=1,W)}.
#' The rest of the terms are estimated using methods specified in \code{method}.
learn_g_tmp <- function(W,
                        A,
                        method,
                        folds,
                        g_bounds,
                        cross_fit_nuisance) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- numeric(length(A))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                         Lrnr_cv_selector$new(loss_loglik_binomial))

    task <- sl3_Task$new(data = data.table(W, A = A), covariates = colnames(W),
                         outcome = "A", outcome_type = "binomial")
    fit <- lrnr$train(task)
    pred <- lrnr$predict(task)

  } else if (method == "glm") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(A[train_idx] ~ ., data = W[train_idx, , drop = FALSE],
                   family = "binomial")
        pred[valid_idx] <<- .bound(as.numeric(predict(fit,
          newdata = W[valid_idx, , drop = FALSE], type = "response"
        )), g_bounds)
      })
    } else {
      # no cross fit
      fit <- glm(A ~ ., data = W, family = "binomial")
      pred <- .bound(as.numeric(predict(fit, newdata = W, type = "response")),
                     g_bounds)
    }

  } else if (method == "glmnet") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(x = W[train_idx, , drop = FALSE], y = A[train_idx],
                         keep = TRUE, alpha = 1, nfolds = length(folds),
                         family = "binomial")
        pred[valid_idx] <<- .bound(as.numeric(predict(fit,
          newx = W[valid_idx, , drop = FALSE], s = "lambda.min", type = "response"
        )), g_bounds)
      })
    } else {
      # no cross fit
      fit <- cv.glmnet(x = W, y = A, keep = TRUE, alpha = 1,
                       nfolds = length(folds), family = "binomial")
      pred <- .bound(as.numeric(predict(fit, newx = W, s = "lambda.min",
                                        type = "response")), g_bounds)
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', 'sl3', or a
         list of sl3 learners.")
  }

  return(pred)
}
