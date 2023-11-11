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
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
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
learn_g <- function(S,
                    W,
                    A,
                    g_rct,
                    controls_only,
                    method,
                    v_folds,
                    g_bounds) {

  if (method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- numeric(length(A))

  if (is.list(method)) {
    # P(S=1|W)
    lrnr_stack <- Stack$new(method)
    lrnr_s_w <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                             Lrnr_cv_selector$new(loss_loglik_binomial))
    task_s_w <- sl3_Task$new(data = data.table(W, S = S),
                             covariates = colnames(W),
                             outcome = "S", outcome_type = "binomial")
    fit_s_w <- lrnr_s_w$train(task_s_w)
    pred_s_w <- fit_s_w$predict(task_s_w)

    if (controls_only) {
      pred <- g_rct*pred_s_w
    } else {
      # P(A=1|S=0,W)
      lrnr_a_ws0 <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                                 Lrnr_cv_selector$new(loss_loglik_binomial))
      task_a_ws0 <- sl3_Task$new(data = data.table(W, A = A)[S == 0,],
                                 covariates = colnames(W),
                                 outcome = "A", outcome_type = "binomial")
      task_a_ws0_pred <- sl3_Task$new(data = data.table(W, A = A),
                                      covariates = colnames(W),
                                      outcome = "A", outcome_type = "binomial")
      fit_a_ws0 <- lrnr_a_ws0$train(task_a_ws0)
      pred_a_ws0 <- fit_a_ws0$predict(task_a_ws0_pred)
      pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
    }

  } else if (method == "glm") {
    # P(S=1|W)
    fit_s_w <- glm(S ~., data = data.frame(W), family = "binomial")
    pred_s_w <- as.numeric(predict(fit_s_w, newdata = data.frame(W),
                                   type = "response"))

    if (controls_only) {
      pred <- g_rct*pred_s_w
    } else {
      # P(A=1|S=0,W)
      fit_a_ws0 <- glm(A[S == 0] ~., data = data.frame(W[S == 0,]),
                       family = "binomial")
      pred_a_ws0 <- as.numeric(predict(fit_a_ws0, newdata = data.frame(W),
                                       type = "response"))
      pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
    }

  } else if (method == "glmnet") {
    # P(S=1|W)
    fit_s_w <- cv.glmnet(x = as.matrix(W), y = S,
                         family = "binomial", keep = TRUE, alpha = 1,
                         nfolds = v_folds)
    pred_s_w <- as.numeric(predict(fit_s_w, newx = as.matrix(W),
                                   s = "lambda.min", type = "response"))

    if (controls_only) {
      pred <- g_rct*pred_s_w
    } else {
      # P(A=1|S=0,W)
      fit_a_ws0 <- cv.glmnet(x = as.matrix(W[S == 0,]), y = A[S == 0],
                             family = "binomial", keep = TRUE, alpha = 1,
                             nfolds = v_folds)
      pred_a_ws0 <- as.numeric(predict(fit_a_ws0, newx = as.matrix(W),
                                       s = "lambda.min", type = "response"))
      pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(.bound(pred, g_bounds))
}
