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
#' W1 <- rnorm(n); W2 <- rnorm(n); W <- cbind(W1, W2)
#' A <- numeric(n)
#' A[S == 1] <- rbinom(sum(S), 1, 0.67)
#' A[S == 0] <- rbinom(n-sum(S), 1, plogis(1.2*W1-0.9*W2))
#'
#' # learn treatment mechanism
#' g <- learn_g(S, W, A, 0.67, FALSE, "glm", 5, c(0.01, 0.99))
learn_g <- function(S,
                    W,
                    A,
                    g_rct,
                    controls_only,
                    method,
                    folds,
                    g_bounds,
                    cross_fit_nuisance) {

  if (is.character(method) && method == "sl3") {
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

    if (controls_only) {
      # control
      X_W <- data.frame(W)

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit_S_W <- glm(S[train_idx] ~., data = X_W[train_idx, ],
                         family = "binomial")
          pred_S_W <- as.numeric(predict(
            fit_S_W, newdata = X_W[valid_idx,], type = "response"))
          pred[valid_idx] <<- g_rct*pred_S_W
        })

      } else {
        # no cross fit
        fit <- glm(S ~., data = X_W, family = "binomial")
        pred_S_W <- as.numeric(predict(fit, newdata = X_W, type = "response"))
        pred <- g_rct*pred_S_W
      }

    } else {
      # control + treatment
      X_W <- data.frame(W)
      X_SW <- as.data.frame(model.matrix(as.formula("~-1+.+S:."),
                                         data = data.frame(S = S, W)))
      X_0W <- as.data.frame(model.matrix(as.formula("~-1+.+S:."),
                                         data = data.frame(S = 0, W)))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit_S_W <- glm(S[train_idx] ~., data = X_W[train_idx,],
                         family = "binomial")
          fit_A_SW <- glm(A[train_idx] ~., data = X_SW[train_idx,],
                          family = "binomial")
          pred_S_W <- as.numeric(predict(
            fit_S_W, newdata = X_W[valid_idx,], type = "response"))
          pred_A_0W <- as.numeric(predict(
            fit_A_SW, newdata = X_0W[valid_idx,], type = "response"))
          pred[valid_idx] <<- g_rct*pred_S_W+pred_A_0W*(1-pred_S_W)
        })

      } else {
        # no cross fit
        fit_S_W <- glm(S ~., data = X_W, family = "binomial")
        pred_S_W <- as.numeric(predict(fit_S_W, newdata = X_W,
                                       type = "response"))
        fit_A_SW <- glm(A ~., data = X_SW, family = "binomial")
        pred_A_0W <- as.numeric(predict(fit_A_SW, newdata = X_0W,
                                        type = "response"))
        pred <- g_rct*pred_S_W+pred_A_0W*(1-pred_S_W)

      }
    }

  } else if (method == "glmnet") {

    if (controls_only) {
      # control
      X_W <- as.matrix(W)

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit_S_W <- cv.glmnet(x = X_W[train_idx,], y = S[train_idx],
                               family = "binomial", keep = TRUE, alpha = 1,
                               nfolds = length(folds))
          pred_S_W <- as.numeric(predict(
            fit_S_W, newx = X_W[valid_idx,], s = "lambda.min",
            type = "response"))
          pred[valid_idx] <<- g_rct*pred_S_W
        })

      } else {
        # no cross fit
        fit_S_W <- cv.glmnet(x = X_W, y = S, family = "binomial", keep = TRUE,
                             alpha = 1, nfolds = length(folds))
        pred_S_W <- as.numeric(predict(
          fit_S_W, newx = X_W, s = "lambda.min", type = "response"))
        pred <- g_rct*pred_S_W
      }

    } else {
      # control + treatment
      X_W <- as.matrix(W)
      X_SW <- model.matrix(as.formula("~-1+.+S:."), data = data.frame(S = S, W))
      X_0W <- model.matrix(as.formula("~-1+.+S:."), data = data.frame(S = 0, W))

      if (cross_fit_nuisance) {
        # cross fit
        walk(folds, function(.x) {
          train_idx <- .x$training_set
          valid_idx <- .x$validation_set
          fit_S_W <- cv.glmnet(x = X_W[train_idx,], y = S[train_idx],
                               family = "binomial", keep = TRUE, alpha = 1,
                               nfolds = length(folds))
          fit_A_SW <- cv.glmnet(x = X_SW[train_idx,], y = A[train_idx],
                                family = "binomial", keep = TRUE, alpha = 1,
                                nfolds = length(folds))
          pred_S_W <- as.numeric(predict(
            fit_S_W, newx = X_W[valid_idx,], s = "lambda.min",
            type = "response"))
          pred_A_0W <- as.numeric(predict(
            fit_A_SW, newx = X_0W[valid_idx,], s = "lambda.min",
            type = "response"))
          pred[valid_idx] <<- g_rct*pred_S_W+pred_A_0W*(1-pred_S_W)
        })

      } else {
        # no cross fit
        fit_S_W <- cv.glmnet(x = X_W, y = S,
                             family = "binomial", keep = TRUE, alpha = 1,
                             nfolds = length(folds))
        fit_A_SW <- cv.glmnet(x = X_SW, y = A,
                              family = "binomial", keep = TRUE, alpha = 1,
                              nfolds = length(folds))
        pred_S_W <- as.numeric(predict(fit_S_W, newx = X_W,
                                       s = "lambda.min", type = "response"))
        pred_A_0W <- as.numeric(predict(fit_A_SW, newx = X_0W,
                                        s = "lambda.min", type = "response"))
        pred <- g_rct*pred_S_W+pred_A_0W*(1-pred_S_W)
      }
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(.bound(pred, g_bounds))
}
