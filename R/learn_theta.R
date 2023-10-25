#' @description Function to learn the conditional mean of outcome given treatment and baseline covariates
#'
#' @param W A matrix of baseline covariates
#' @param A A vector of treatment indicators, A = 1 for treatment, A = 0 for control
#' @param Y A vector of outcomes
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_theta <- function(W, A, Y,
                        controls_only,
                        method,
                        v_folds) {
  if (method == "sl3") {
    method <- get_default_sl3_learners("gaussian")
  }

  pred <- numeric(length = length(A))

  if (method == "glm") {
    # A = 0
    fit_A0 <- glm(Y[A == 0] ~., data = data.frame(W[A == 0, ]), family = "gaussian")
    A0 <- as.numeric(predict(fit_A0, newdata = data.frame(W[A == 0, ]), type = "response"))
    pred[A == 0] <- A0

    if (!controls_only) {
      # A = 1
      fit_A1 <- glm(Y[A == 1] ~., data = data.frame(W[A == 1, ]), family = "gaussian")
      A1 <- as.numeric(predict(fit_A1, newdata = data.frame(W[A == 1, ]), type = "response"))
      pred[A == 1] <- A1
    }

  } else if (method == "glmnet") {
    # A = 0
    fit_A0 <- cv.glmnet(x = as.matrix(W[A == 0, ]), y = Y[A == 0],
                        keep = TRUE, alpha = 1, nfolds = v_folds, family = "gaussian")
    A0 <- as.numeric(predict(fit, newx = as.matrix(W[A == 0, ]), s = "lambda.min", type = "response"))
    pred[A == 0] <- A0

    if (!controls_only) {
      # A = 1
      fit_A1 <- cv.glmnet(x = as.matrix(W[A == 1, ]), y = Y[A == 1],
                          keep = TRUE, alpha = 1, nfolds = v_folds, family = "gaussian")
      A1 <- as.numeric(predict(fit, newx = as.matrix(W[A == 1, ]), s = "lambda.min", type = "response"))
      pred[A == 1] <- A1
    }

  } else if (is.list(method)) {
    lrnr_stack <- Stack$new(method)
    lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack),
                               Lrnr_cv_selector$new(loss_squared_error))

    if (controls_only) {
      task_theta <- sl3_Task$new(data = data.table(W, Y = Y)[A == 0, ],
                                 covariates = colnames(W),
                                 outcome = "Y", outcome_type = "continuous")
      fit_theta <- lrnr_theta$train(task_theta)
      pred[A == 0] <- fit_theta$predict(task_theta)

    } else {
      task_theta <- sl3_Task$new(data = data.table(W, Y = Y, A = A),
                                 covariates = c(colnames(W), "A"),
                                 outcome = "Y", outcome_type = "continuous")
      fit_theta <- lrnr_theta$train(task_theta)
      pred <- fit_theta$predict(task_theta)
    }
  }

  return(list(pred = pred))
}
