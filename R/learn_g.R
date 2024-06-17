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
#' @param data A `data.table` object containing the data.
#' @param W A `character` vector of baseline covariates.
#' @param A A `character` string of the treatment variable.
#' @param method Learning method. `glm` for main-term linear model,
#' `glmnet` for main-term lasso, `sl3` for super learner-based estimation
#' using default learners, or a `list` of `sl3` learners for super
#' learner-based estimation with custom learners.
#' @param folds A `list` of cross-validation folds made using `origami` package
#' (to be used whenever necessary).
#' @param g_bounds A numeric vector of lower and upper bounds for the
#' treatment mechanism. The first element is the lower bound, and the second
#' element is the upper bound.
#'
#' @returns A numeric vector of estimated treatment probabilities.
learn_g <- function(data,
                    W,
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
    lrnrs <- make_learner(Pipeline,
                          Lrnr_cv$new(lrnr_stack),
                          Lrnr_cv_selector$new(loss_loglik_binomial))
    task <- sl3_Task$new(data = data,
                         covariates = W,
                         outcome = A,
                         outcome_type = "binomial")
    fit <- lrnrs$train(task)
    pred <- fit$predict(task)

  } else if (method == "glm") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(data[[A]][train_idx] ~ .,
                   data = data[train_idx, ..W],
                   family = "binomial")
        pred[valid_idx] <<- predict(fit, newdata = data[valid_idx, ..W],
                                    type = "response")
      })
    } else {
      # no cross fit
      fit <- glm(data[[A]] ~ ., data = data[, ..W], family = "binomial")
      pred <- predict(fit, newdata = data[, ..W], type = "response")
    }

  } else if (method == "glmnet") {

    if (cross_fit_nuisance) {
      # cross fit
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(x = as.matrix(data[train_idx, ..W]),
                         y = data[[A]][train_idx],
                         family = "binomial", keep = TRUE, alpha = 1,
                         nfolds = length(folds))
        pred[valid_idx] <<- predict(fit,
                                    newx = as.matrix(data[valid_idx, ..W]),
                                    s = "lambda.min", type = "response")
      })
    } else {
      # no cross fit
      fit <- cv.glmnet(x = data.table:::as.matrix(data[, ..W]), y = data[[A]],
                       family = "binomial", keep = TRUE, alpha = 1,
                       nfolds = length(folds))
      pred <- predict(fit, newx = data.table:::as.matrix(data[, ..W]),
                      s = "lambda.min", type = "response")
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(.bound(as.numeric(pred), g_bounds))
}
