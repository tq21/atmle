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
#' @returns A numeric vector of the estimated values.
#'
#' @examples
#' # simulate data
#' set.seed(123)
#' n <- 500
#' S <- rbinom(n, 1, 0.5)
#' W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
#' A <- numeric(n)
#' A[S == 1] <- rbinom(sum(S), 1, 0.67)
#' A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
#' UY <- rnorm(n, 0, 1)
#' U_bias <- rnorm(n, 0, 0.5)
#' Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)+(1-S)*U_bias
#'
#' theta <- learn_theta(W, A, Y, FALSE, "glm", 5, "gaussian", c(-4, 4))
learn_theta <- function(W,
                        A,
                        Y,
                        controls_only,
                        method,
                        folds,
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

  pred <- rep(NA, length(A))

  if (is.list(method)) {
    if (family == "gaussian") {
      lrnr_stack <- Stack$new(method)
      lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                 Lrnr_cv_selector$new(loss_squared_error))

      if (controls_only) {
        task_theta <- sl3_Task$new(data = data.table(W, Y = Y)[A == 0, ],
                                   covariates = colnames(W),
                                   outcome = "Y", outcome_type = "continuous")
        fit_theta <- lrnr_theta$train(task_theta)
        pred[A == 0] <- .bound(fit_theta$predict(task_theta), theta_bounds)

      } else {
        task_theta <- sl3_Task$new(data = data.table(W, Y = Y, A = A),
                                   covariates = c(colnames(W), "A"),
                                   outcome = "Y", outcome_type = "continuous")
        fit_theta <- lrnr_theta$train(task_theta)
        pred <- .bound(fit_theta$predict(task_theta), theta_bounds)
      }
    } else if (family == "binomial") {
      lrnr_stack <- Stack$new(method)
      lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                 Lrnr_cv_selector$new(loss_loglik_binomial))

      if (controls_only) {
        task_theta <- sl3_Task$new(data = data.table(W, Y = Y)[A == 0, ],
                                   covariates = colnames(W),
                                   outcome = "Y", outcome_type = "binomial")
        fit_theta <- lrnr_theta$train(task_theta)
        pred[A == 0] <- .bound(fit_theta$predict(task_theta), theta_bounds)

      } else {
        task_theta <- sl3_Task$new(data = data.table(W, Y = Y, A = A),
                                   covariates = c(colnames(W), "A"),
                                   outcome = "Y", outcome_type = "binomial")
        fit_theta <- lrnr_theta$train(task_theta)
        pred <- .bound(fit_theta$predict(task_theta), theta_bounds)
      }
    } else {
      stop("Invalid family. Must be either 'gaussian' or 'binomial'.")
    }

  } else if (method == "glm") {
    if (controls_only) {
      fit <- glm(Y[A == 0] ~., data = data.frame(W[A == 0, ]), family = family)
      A0 <- .bound(as.numeric(predict(fit, newdata = data.frame(W[A == 0, ]),
                                      type = "response")), theta_bounds)
      pred[A == 0] <- A0

    } else {
      X <- as.data.frame(model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A)))
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- glm(Y[train_idx] ~., data = X[train_idx, ], family = family)
        pred[valid_idx] <<- .bound(as.numeric(predict(fit, newdata = X[valid_idx,],
                                                      type = "response")), theta_bounds)
      })
    }

  } else if (method == "glmnet") {
    if (controls_only) {
      fit_A0 <- cv.glmnet(x = as.matrix(W[A == 0, ]), y = Y[A == 0],
                          keep = TRUE, alpha = 1, nfolds = length(folds),
                          family = family)
      A0 <- .bound(as.numeric(predict(fit_A0, newx = as.matrix(W[A == 0, ]),
                                      s = "lambda.min", type = "response")),
                   theta_bounds)
      pred[A == 0] <- A0

    } else {
      X <- as.data.frame(model.matrix(as.formula("~-1+.+A:."), data = data.frame(W, A = A)))
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set
        fit <- cv.glmnet(x = as.matrix(X[train_idx, ]), y = Y[train_idx],
                         keep = TRUE, alpha = 1, nfolds = length(folds), family = family)
        pred[valid_idx] <<- .bound(as.numeric(predict(fit, newx = as.matrix(X[valid_idx, ]),
                                                    s = "lambda.min", type = "response")),
                                   theta_bounds)
      })
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(pred)
}
