#' Learn nuisance function: conditional mean of outcome given baseline
#' covariates
#'
#' @description Function to learn the conditional mean of outcome given
#' baseline covariates, \eqn{\tilde{\theta}(W)=\mathbb{E}(Y\mid W)}.
#'
#' @export
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
#'
#' @param W A matrix of baseline covariates.
#' @param Y A vector of outcomes.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#' @param family A character string specifying the family of the outcome
#' \eqn{Y}.
#'
#' @returns A numeric vector of the estimated values.
learn_theta_tilde <- function(W, Y,
                              method,
                              v_folds,
                              family,
                              theta_bound) {

  if (method == "sl3") {
    method <- get_default_sl3_learners()
  }

  if (is.null(theta_bound)) {
    if (family == "gaussian") {
      theta_bound <- c(min(Y), max(Y))
    } else if (family == "binomial") {
      theta_bound <- c(0.01, 0.99)
    }
  }

  pred <- numeric(length(Y))

  if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_theta_tilde <- NULL
    task_theta_tilde <- NULL
    if (family == "gaussian") {
      lrnr_theta_tilde <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                       Lrnr_cv_selector$new(loss_squared_error))
      task_theta_tilde <- sl3_Task$new(data = data.table(W, Y = Y),
                                       covariates = colnames(W),
                                       outcome = "Y", outcome_type = "continuous")
    } else if (family == "binomial") {
      lrnr_theta_tilde <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                       Lrnr_cv_selector$new(loss_loglik_binomial))
      task_theta_tilde <- sl3_Task$new(data = data.table(W, Y = Y),
                                       covariates = colnames(W),
                                       outcome = "Y", outcome_type = "binomial")
    }

    fit_theta_tilde <- lrnr_theta_tilde$train(task_theta_tilde)
    pred <- fit_theta_tilde$predict(task_theta_tilde)
  } else if (method == "glm") {
    fit <- glm(Y ~., data = data.frame(W), family = family)
    pred <- as.numeric(predict(fit, newdata = data.frame(W), type = "response"))

  } else if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W), y = Y,
                     keep = TRUE, alpha = 1, nfolds = v_folds, family = family)
    pred <- as.numeric(predict(fit, newx = as.matrix(W), s = "lambda.min", type = "response"))
  }

  print("low theta_tilde: " %+% min(pred))
  print("high theta_tilde: " %+% max(pred))

  return(.bound(pred, theta_bound))
}
