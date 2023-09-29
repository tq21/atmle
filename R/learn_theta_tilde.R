#' @description Function to learn the conditional mean of outcome given baseline covariates
#'
#' @param W A matrix of baseline covariates
#' @param Y A vector of outcomes
#' @param family "gaussian" or "binomial"
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_theta_tilde <- function(W, Y, family = "gaussian", method = "glmnet", v_folds = 5) {
  pred <- NULL

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W), y = Y,
                     keep = TRUE, alpha = 1, nfolds = v_folds, family = family)
    pred <- as.numeric(predict(fit, newx = as.matrix(W), s = "lambda.min", type = "response"))
  } else if (method == "glm") {
    fit <- glm(Y ~., data = data.frame(W), family = family)
    pred <- as.numeric(predict(fit, newdata = data.frame(W), type = "response"))
  } else if (method == "sl3") {
    lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 3, family = "binomial"),
                                 Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0)))
    lrnr_theta_tilde <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new(loss_squared_error))
    task_theta_tilde <- sl3_Task$new(data = data.table(W, Y = Y),
                                     covariates = colnames(W),
                                     outcome = "Y", outcome_type = "continuous")
    fit_theta_tilde <- lrnr_theta_tilde$train(task_theta_tilde)
    pred <- fit_theta_tilde$predict(task_theta_tilde)
  }

  return(pred)
}
