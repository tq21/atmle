#' @description Function to learn the treatment mechanism
#'
#' @param S A vector of study indicators, S = 1 for RCT, S = 0 for RWD
#' @param W A matrix of baseline covariates
#' @param A A vector of treatment indicators, A = 1 for treatment, A = 0 for control
#' @param g_rct A vector of treatment probability in RCT
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_g <- function(S, W, A,
                    g_rct,
                    controls_only,
                    method,
                    v_folds = 5) {

  pred <- numeric(length(A))

  if (method == "glm") {
    # P(S=1|W)
    fit_s_w <- glm(S ~., data = data.frame(W), family = "binomial")
    pred_s_w <- as.numeric(predict(fit_s_w, newdata = data.frame(W), type = "response"))

    if (controls_only) {
      pred <- g_rct*pred_s_w
    } else {
      # P(A=1|S=0,W)
      fit_a_ws0 <- glm(A[S == 0] ~., data = data.frame(W[S == 0,]), family = "binomial")
      pred_a_ws0 <- as.numeric(predict(fit_a_ws0, newdata = data.frame(W), type = "response"))
      pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
    }

  } else if (method == "glmnet") {
    # P(S=1|W)
    fit_s_w <- cv.glmnet(x = as.matrix(W), y = S,
                         family = "binomial", keep = TRUE, alpha = 1, nfolds = v_folds)
    pred_s_w <- as.numeric(predict(fit_s_w, newx = as.matrix(W), s = "lambda.min", type = "response"))

    if (controls_only) {
      pred <- g_rct*pred_s_w
    } else {
      # P(A=1|S=0,W)
      fit_a_ws0 <- cv.glmnet(x = as.matrix(W[S == 0,]), y = A[S == 0],
                             family = "binomial", keep = TRUE, alpha = 1, nfolds = v_folds)
      pred_a_ws0 <- as.numeric(predict(fit_a_ws0, newx = as.matrix(W), s = "lambda.min", type = "response"))
      pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
    }

  } else if (is.list(method)) {
    # P(S=1|W)
    lrnr_stack <- Stack$new(method)
    lrnr_s_w <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                             Lrnr_cv_selector$new(loss_loglik_binomial))
    task_s_w <- sl3_Task$new(data = data.table(W, S = S),
                             covariates = colnames(W),
                             outcome = "S", outcome_type = "binomial")
    fit_s_w <- lrnr_s$train(task_s_w)
    pred_s_w <- fit_s$predict(task_s_w)

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
  }

  return(pred)
}
