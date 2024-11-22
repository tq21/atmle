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
learn_g <- function(S,
                    W,
                    A,
                    method,
                    controls_only,
                    folds,
                    g_bounds,
                    cross_fit_nuisance) {

  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred_S_W <- pred_A_S1_W <- pred_A_S0_W <- numeric(length(S))

  if (is.list(method)) {
    # TODO: factorization
    lrnr_stack <- Stack$new(method)
    lrnr <- make_learner(
      Pipeline, Lrnr_cv$new(lrnr_stack),
      Lrnr_cv_selector$new(loss_loglik_binomial)
    )

    # P(S=1|W)
    task_S_W <- sl3_Task$new(
      data = data.table(W, S = S),
      covariates = colnames(W),
      outcome = "S", outcome_type = "binomial"
    )
    fit_S_W <- lrnr$train(task_S_W)
    pred_S_W <- fit_S_W$predict(task_S_W)

    # P(A=1|S=1,W)
    task_A_S1_W <- sl3_Task$new(
      data = data.table(W[S == 1, ], A = A[S == 1]),
      covariates = colnames(W),
      outcome = "A", outcome_type = "binomial"
    )
    fit_A_S1_W <- lrnr$train(task_A_S1_W)
    pred_A_S1_W <- fit_A_S1_W$predict(task_S_W)

    if (!controls_only) {
      # P(A=1|S=0,W)
      task_A_S0_W <- sl3_Task$new(
        data = data.table(W[S == 0, ], A = A[S == 0]),
        covariates = colnames(W),
        outcome = "A", outcome_type = "binomial"
      )
      fit_A_S0_W <- lrnr$train(task_A_S0_W)
      pred_A_S0_W <- fit_A_S0_W$predict(task_S_W)
    }

  } else if (method == "glm") {

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set

        # P(S=1|W)
        fit_S_W <- glm(S[train_idx] ~ ., data = W[train_idx, ], family = "binomial")
        pred_S_W[valid_idx] <<- as.numeric(predict(fit_S_W, newdata = W[valid_idx, ], type = "response"))

        # P(A=1|S=1,W)
        fit_A_S1_W <- glm(A[train_idx][S[train_idx] == 1] ~ .,
                          data = W[train_idx, ][S[train_idx] == 1, ], family = "binomial")
        pred_A_S1_W[valid_idx] <<- as.numeric(predict(fit_A_S1_W, newdata = W[valid_idx, ], type = "response"))

        if (!controls_only) {
          # P(A=1|S=0,W)
          fit_A_S0_W <- glm(A[train_idx][S[train_idx] == 0] ~ .,
                            data = W[train_idx, ][S[train_idx] == 0, ], family = "binomial")
          pred_A_S0_W[valid_idx] <<- as.numeric(predict(fit_A_S0_W, newdata = W[valid_idx, ], type = "response"))
        }
      })

    } else {
      # P(S=1|W)
      fit_S_W <- glm(S ~ ., data = W, family = "binomial")
      pred_S_W <- as.numeric(predict(fit_S_W, newdata = W, type = "response"))

      # P(A=1|S=1,W)
      fit_A_S1_W <- glm(A[S == 1] ~ ., data = W[S == 1, ], family = "binomial")
      pred_A_S1_W <- as.numeric(predict(fit_A_S1_W, newdata = W, type = "response"))

      if (!controls_only) {
        # P(A=1|S=0,W)
        fit_A_S0_W <- glm(A[S == 0] ~ ., data = W[S == 0, ], family = "binomial")
        pred_A_S0_W <- as.numeric(predict(fit_A_S0_W, newdata = W, type = "response"))
      }
    }

  } else if (method == "glmnet") {

    if (cross_fit_nuisance) {
      walk(folds, function(.x) {
        train_idx <- .x$training_set
        valid_idx <- .x$validation_set

        # P(S=1|W)
        fit_S_W <- cv.glmnet(x = as.matrix(W[train_idx, , drop = FALSE]),
                             y = S[train_idx], family = "binomial", keep = TRUE, alpha = 1,
                             nfolds = length(folds))
        pred_S_W[valid_idx] <<- .bound(as.numeric(predict(
          fit_S_W,
          newx = as.matrix(W[valid_idx, , drop = FALSE]), s = "lambda.min",
          type = "response"
        )), g_bounds)

        # P(A=1|S=1,W)
        fit_A_S1_W <- cv.glmnet(x = as.matrix(W[train_idx, ][S[train_idx] == 1, ]),
                                y = A[train_idx][S[train_idx] == 1],
                                family = "binomial", keep = TRUE, alpha = 1,
                                nfolds = length(folds))
        pred_A_S1_W[valid_idx] <<- as.numeric(predict(fit_A_S1_W,
                                                     newx = as.matrix(W[valid_idx, ]),
                                                     s = "lambda.min",
                                                     type = "response"))

        if (!controls_only) {
          # P(A=1|S=0,W)
          fit_A_S0_W <- cv.glmnet(x = as.matrix(W[train_idx, ][S[train_idx] == 0, ]),
                                  y = A[train_idx][S[train_idx] == 0],
                                  family = "binomial", keep = TRUE, alpha = 1,
                                  nfolds = length(folds))
          pred_A_S0_W[valid_idx] <<- as.numeric(predict(fit_A_S0_W,
                                                       newx = as.matrix(W[valid_idx, ]),
                                                       s = "lambda.min",
                                                       type = "response"))
        }
      })
    } else {

      # P(S=1|W)
      fit_S_W <- cv.glmnet(x = as.matrix(W), y = S, family = "binomial",
                           keep = TRUE, alpha = 1, nfolds = v_folds)
      pred_S_W <- as.numeric(predict(fit_S_W, newx = as.matrix(W), s = "lambda.min", type = "response"))

      # P(A=1|S=1,W)
      fit_A_S1_W <- cv.glmnet(x = as.matrix(W[S == 1, ]), y = A[S == 1], family = "binomial",
                              keep = TRUE, alpha = 1, nfolds = v_folds)
      pred_A_S1_W <- as.numeric(predict(fit_A_S1_W, newx = as.matrix(W), s = "lambda.min", type = "response"))

      if (!controls_only) {
        # P(A=1|S=0,W)
        fit_A_S0_W <- cv.glmnet(x = as.matrix(W[S == 0, ]), y = A[S == 0], family = "binomial",
                                keep = TRUE, alpha = 1, nfolds = v_folds)
        pred_A_S0_W <- as.numeric(predict(fit_A_S0_W, newx = as.matrix(W), s = "lambda.min", type = "response"))
      }
    }

  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', 'sl3', or a
         list of sl3 learners.")
  }

  # P(A=1|W)=P(A=1|S=1,W)P(S=1|W)+P(A=1|S=0,W)P(S=0|W)
  pred <- pred_A_S1_W * pred_S_W + pred_A_S0_W * (1 - pred_S_W)

  return(list(pred = .bound(pred, g_bounds),
              pred_S_W = pred_S_W,
              pred_A_S1_W = pred_A_S1_W,
              pred_A_S0_W = pred_A_S0_W))
}
