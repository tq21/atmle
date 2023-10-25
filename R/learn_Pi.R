#' @description Function to learn the trial enrollment mechanism
#'
#' @param S A vector of study indicators, S = 1 for RCT, S = 0 for RWD
#' @param W A matrix of baseline covariates
#' @param A A vector of treatment indicators, A = 1 for treatment, A = 0 for control
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_Pi <- function(S, W, A,
                     controls_only,
                     method,
                     v_folds) {
  if (method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- numeric(length(A))
  A1 <- numeric(length(A))
  A0 <- numeric(length(A))

  if (method == "glm") {
    # A = 0
    fit_A0 <- glm(S[A == 0] ~., data = data.frame(W[A == 0,]), family = "binomial")
    A0 <- as.numeric(predict(fit_A0, newdata = data.frame(W), type = "response"))
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1 # no treated in external
    } else {
      # A = 1
      fit_A1 <- glm(S[A == 1] ~., data = data.frame(W[A == 1,]), family = "binomial")
      A1 <- as.numeric(predict(fit_A1, newdata = data.frame(W), type = "response"))
      pred[A == 1] <- A1[A == 1]
    }

  } else if (method == "glmnet") {
    # A = 0
    fit_A0 <- cv.glmnet(x = as.matrix(W[A == 0, ]), y = S[A == 0],
                        keep = TRUE, alpha = 1, nfolds = v_folds, family = "binomial")
    A0 <- as.numeric(predict(fit_A0, newx = as.matrix(W), s = "lambda.min", type = "response"))
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1
    } else {
      # A = 1
      fit_A1 <- cv.glmnet(x = as.matrix(W[A == 1, ]), y = S[A == 1],
                          keep = TRUE, alpha = 1, nfolds = v_folds, family = "binomial")
      A1 <- as.numeric(predict(fit_A1, newx = as.matrix(W), s = "lambda.min", type = "response"))
      pred[A == 1] <- A1[A == 1]
    }

  } else if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_Pi <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                            Lrnr_cv_selector$new(loss_loglik_binomial))
    task_Pi <- sl3_Task$new(data = data.table(W, A = A, S = S),
                            covariates = c(colnames(W), "A"),
                            outcome = "S", outcome_type = "binomial")
    fit_Pi <- lrnr_Pi$train(task_Pi)
    task_Pi_A0 <- sl3_Task$new(data = data.table(W, A = 0, S = S),
                               covariates = c(colnames(W), "A"),
                               outcome = "S", outcome_type = "binomial")
    A0 <- fit_Pi$predict(task_Pi_A0)
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1 # no treated in external
    } else {
      task_Pi_A1 <- sl3_Task$new(data = data.table(W, A = 1, S = S),
                                 covariates = c(colnames(W), "A"),
                                 outcome = "S", outcome_type = "binomial")
      A1 <- fit_Pi$predict(task_Pi_A1)
      pred[A == 1] <- A1[A == 1]
    }
  }

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}
