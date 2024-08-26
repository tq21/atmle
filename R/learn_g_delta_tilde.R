#' @title Learn nuisance function: conditional probability of observing the
#' outcome given baseline covariates and treatment
#'
#' @description Function to learn the conditional probability of observing the
#' outcome given baseline covariates and treatment,
#' \eqn{\tilde{g}_\Delta(1\mid W,A)=\mathbb{P}(\Delta=1\mid W,A)}. Only applicable when
#' the outcome is subject to missingness.
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
#' @param delta A vector of missingness indicators.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param g_bounds A numeric vector of lower and upper bounds for the
#' missingness mechanism. The first element is the lower bound, and the second
#' element is the upper bound.
#'
#' @returns A numeric vector of the estimated values.
learn_g_delta_tilde <- function(W,
                                A,
                                delta,
                                method,
                                folds,
                                g_bounds,
                                cross_fit_nuisance = FALSE) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- A0 <- A1 <- numeric(length(A))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_delta <- make_learner(
      Pipeline, Lrnr_cv$new(lrnr_stack),
      Lrnr_cv_selector$new(loss_loglik_binomial)
    )
    task_train <- sl3_Task$new(
      data = data.table(W, delta = delta, A = A),
      covariates = c(colnames(W), "A"),
      outcome = "delta", outcome_type = "binomial"
    )
    task_A0 <- sl3_Task$new(
      data = data.table(W, delta = delta, A = 0),
      covariates = c(colnames(W), "A"),
      outcome = "delta", outcome_type = "binomial"
    )
    task_A1 <- sl3_Task$new(
      data = data.table(W, delta = delta, A = 1),
      covariates = c(colnames(W), "A"),
      outcome = "delta", outcome_type = "binomial"
    )
    fit_delta <- lrnr_delta$train(task_train)
    pred <- fit_delta$predict(task_train)
    A0 <- fit_delta$predict(task_A0)
    A1 <- fit_delta$predict(task_A1)
  } else if (method == "glm") {
    X <- as.data.frame(model.matrix(
      as.formula("~-1+.+A:."),
      data = data.frame(W, A = A)
    ))
    X_A0 <- as.data.frame(model.matrix(
      as.formula("~-1+.+A:."),
      data = data.frame(W, A = 0)
    ))
    X_A1 <- as.data.frame(model.matrix(
      as.formula("~-1+.+A:."),
      data = data.frame(W, A = 1)
    ))

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(delta[train_idx] ~ ., data = X[train_idx, ], family = "binomial")
        pred[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X[valid_idx, ], type = "response"
        )), g_bounds)
        A0[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X_A0[valid_idx, ], type = "response"
        )), g_bounds)
        A1[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newdata = X_A1[valid_idx, ], type = "response"
        )), g_bounds)
      })
    } else {
      fit <- glm(delta ~ ., data = X, family = "binomial")
      pred <- as.numeric(predict(fit, newdata = X, type = "response"))
      A0 <- as.numeric(predict(fit, newdata = X_A0, type = "response"))
      A1 <- as.numeric(predict(fit, newdata = X_A1, type = "response"))
    }
  } else if (method == "glmnet") {
    X <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A))
    X_A0 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 0))
    X_A1 <- model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = 1))

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(
          x = X[train_idx, ],
          y = delta[train_idx],
          keep = TRUE, alpha = 1, nfolds = length(folds),
          family = "binomial"
        )
        pred[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newx = X[valid_idx, ], s = "lambda.min", type = "response"
        )), g_bounds)
        A0[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newx = X_A0[valid_idx, ], s = "lambda.min", type = "response"
        )), g_bounds)
        A1[valid_idx] <<- .bound(as.numeric(predict(
          fit,
          newx = X_A1[valid_idx, ], s = "lambda.min", type = "response"
        )), g_bounds)
      })
    } else {
      fit <- cv.glmnet(
        x = X, y = delta, keep = TRUE, alpha = 1,
        nfolds = length(folds), family = "binomial"
      )
      pred <- as.numeric(predict(
        fit,
        newx = X, s = "lambda.min", type = "response"
      ))
      A0 <- as.numeric(predict(
        fit,
        newx = X_A0, s = "lambda.min", type = "response"
      ))
      A1 <- as.numeric(predict(
        fit,
        newx = X_A1, s = "lambda.min", type = "response"
      ))
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(list(
    pred = .bound(pred, g_bounds),
    A0 = .bound(A0, g_bounds),
    A1 = .bound(A1, g_bounds)
  ))
}
