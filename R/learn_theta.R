#' @description Function to learn the conditional mean of outcome given treatment and baseline covariates
#'
#' @param W A matrix of baseline covariates
#' @param A A vector of treatment indicators, A = 1 for treatment, A = 0 for control
#' @param Y A vector of outcomes
#' @param family "gaussian" or "binomial"
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_theta <- function(W, A, Y, family = "gaussian", method = "glmnet", v_folds = 5) {
  pred <- NULL

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(data.frame(W, A)), y = Y,
                     keep = TRUE, alpha = 1, nfolds = v_folds, family = family)
    pred <- as.numeric(predict(fit, newx = as.matrix(data.frame(W, A)), s = "lambda.min", type = "response"))
  } else if (method == "glm") {
    fit <- glm(Y ~., data = data.frame(W, A), family = family)
    pred <- as.numeric(predict(fit, newdata = data.frame(W, A), type = "response"))
  } else if (method == "sl3") {
    lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 3, family = "binomial"),
                                 Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0)))
    lrnr_theta <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new(loss_squared_error))
    task_theta <- sl3_Task$new(data = data.table(W, A = A, Y = Y),
                               covariates = c(colnames(W), "A"),
                               outcome = "Y", outcome_type = "continuous")
    fit_theta <- lrnr_theta$train(task_theta)
    pred <- fit_theta$predict(task_theta)
  }

  return(list(pred = pred))
}
