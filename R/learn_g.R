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
#'
#' @examples
#' # simulate data
#' set.seed(123)
#' n <- 500
#' S <- rbinom(n, 1, 0.5)
#' W1 <- rnorm(n)
#' W2 <- rnorm(n)
#' W <- cbind(W1, W2)
#' A <- numeric(n)
#' A[S == 1] <- rbinom(sum(S), 1, 0.67)
#' A[S == 0] <- rbinom(n - sum(S), 1, plogis(1.2 * W1 - 0.9 * W2))
#'
#' # learn treatment mechanism
#' g <- learn_g(S, W, A, 0.67, FALSE, "glm", 5, c(0.01, 0.99))
learn_g <- function(W,
                    A,
                    method,
                    folds,
                    g_bounds,
                    cross_fit_nuisance = TRUE) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- numeric(length(A))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr <- make_learner(
      Pipeline, Lrnr_cv$new(lrnr_stack),
      Lrnr_cv_selector$new(loss_loglik_binomial)
    )

    task <- sl3_Task$new(
      data = data.table(W, A = A),
      covariates = colnames(W),
      outcome = "A", outcome_type = "binomial"
    )

    fit <- lrnr$train(task)
    pred <- fit$predict(task)

  } else if (method == "glm") {

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(A[train_idx] ~ ., data = W[train_idx, ], family = "binomial")
        pred[valid_idx] <<- as.numeric(predict(fit, newdata = W[valid_idx, ], type = "response"))
      })
    } else {
      fit <- glm(A ~ ., data = W, family = "binomial")
      pred <- as.numeric(predict(fit, newdata = W, type = "response"))
    }

  } else if (method == "glmnet") {

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fot <- cv.glmnet(x = W[train_idx, ], y = A[train_idx], family = "binomial",
                         keep = TRUE, alpha = 1, nfolds = length(folds))
        pred[valid_idx] <<- as.numeric(predict(fit, newx = W[valid_idx, ], s = "lambda.min", type = "response"))
      })
    } else {
      fit <- cv.glmnet(x = W, y = A, family = "binomial",
                       keep = TRUE, alpha = 1, nfolds = length(folds))
      pred <- as.numeric(predict(fit, newx = W, s = "lambda.min", type = "response"))
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', 'sl3', or a
         list of sl3 learners.")
  }

  return(.bound(pred, g_bounds))
}
