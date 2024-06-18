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
#' @param data A `data.table` object containing the data.
#' @param W A `character` vector of baseline covariate names.
#' @param Y A `character` string of the outcome variable name.
#' @param delta A `numeric` vector of censoring indicators.
#' @param method Learning method. `glm` for main-term linear model,
#' `glmnet` for main-term lasso, `sl3` for super learner-based estimation
#' using default learners, or a `list` of `sl3` learners for super
#' learner-based estimation with custom learners.
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}. Either `gaussian` or `binomial`.
#' @param theta_bounds A numeric vector of lower and upper bounds for the
#' conditional mean of outcome given baseline covariates. The first element is
#' the lower bound, and the second element is the upper bound.
#' @param folds A `list` of cross-validation folds made using `origami` package
#' (to be used whenever necessary).
#' @param cross_fit_nuisance A `logical` indicating whether to use
#' out-of-sample predictions for the nuisance estimates. Default is `TRUE`.
#'
#' @returns A `numeric` vector of the estimated values.
learn_theta_tilde <- function(data,
                              W,
                              Y,
                              delta,
                              method,
                              family,
                              theta_bounds,
                              folds,
                              cross_fit_nuisance = TRUE) {
  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners(family)
  }

  pred <- numeric(length(data[, .N]))

  if (is.null(theta_bounds)) {
    if (family == "gaussian") {
      theta_bounds <- c(-Inf, Inf)
    } else if (family == "binomial") {
      theta_bounds <- c(0.01, 0.99)
    }
  }

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    data_Y_0 <- copy(data); set(data_Y_0, j = Y, value = 0) # this is a trick to mute warning of sl3 task for missing Y
    if (family == "gaussian") {
      lrnr <- make_learner(Pipeline,
                           Lrnr_cv$new(lrnr_stack),
                           Lrnr_cv_selector$new(loss_squared_error))
      task_train <- sl3_Task$new(data = data[delta == 1],
                                 covariates = W,
                                 outcome = Y,
                                 outcome_type = "continuous")
      task_pred <- sl3_Task$new(data = data_Y_0, # dummy outcome so sl3 task does not throw warning
                                covariates = W,
                                outcome = Y,
                                outcome_type = "continuous")
    } else if (family == "binomial") {
      lrnr <- make_learner(Pipeline,
                           Lrnr_cv$new(lrnr_stack),
                           Lrnr_cv_selector$new(loss_loglik_binomial))
      task_train <- sl3_Task$new(data = data[delta == 1],
                                 covariates = W,
                                 outcome = Y,
                                 outcome_type = "binomial")
      task_pred <- sl3_Task$new(data = data_Y_0, # dummy outcome so sl3 task does not throw warning
                                covariates = W,
                                outcome = Y,
                                outcome_type = "binomial")
    }

    fit <- lrnr$train(task_train)
    pred <- fit$predict(task_pred)

  } else if (method == "glm") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(data[train_idx][delta[train_idx] == 1][[Y]] ~ .,
                   data = data[train_idx, ..W][delta[train_idx] == 1],
                   family = family)
        pred[valid_idx] <<- predict(fit, newdata = data[valid_idx, ..W],
                                    type = "response")
      })
    } else {
      # no cross fit
      fit <- glm(data[delta == 1][[Y]] ~ ., data = data[delta == 1, ..W], family = family)
      pred <- predict(fit, newdata = data[, ..W], type = "response")
    }

  } else if (method == "glmnet") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(x = as.matrix(data[train_idx, ..W][delta[train_idx] == 1]),
                         y = data[train_idx][delta[train_idx] == 1][[Y]],
                         keep = TRUE, alpha = 1, nfolds = length(folds), family = family)
        pred[valid_idx] <<- predict(fit, newx = as.matrix(data[valid_idx, ..W]),
                                    s = "lambda.min", type = "response")
      })
    } else {
      # no cross fit
      fit <- cv.glmnet(x = as.matrix(data[delta == 1, ..W]),
                       y = data[delta == 1][[Y]], keep = TRUE, alpha = 1,
                       nfolds = length(folds), family = family)
      pred <- predict(fit, newx = data[, ..W], s = "lambda.min", type = "response")
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(.bound(as.numeric(pred), theta_bounds))
}
