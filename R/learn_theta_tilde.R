#' @description Function to learn the conditional mean of outcome given baseline covariates
#'
#' @param W A matrix of baseline covariates
#' @param Y A vector of outcomes
#' @param family "gaussian" or "binomial"
#' @param method Learning method,
#'               "glm" for linear regression,
#'               "glmnet" for lasso,
#'               "sl3" for super learner
learn_theta_tilde <- function(W, Y,
                              method,
                              v_folds,
                              family) {

  if (method == "sl3") {
    method <- get_default_sl3_learners("gaussian")
  }

  pred <- numeric(length(Y))

  if (method == "glm") {
    fit <- glm(Y ~., data = data.frame(W), family = family)
    pred <- as.numeric(predict(fit, newdata = data.frame(W), type = "response"))

  } else if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W), y = Y,
                     keep = TRUE, alpha = 1, nfolds = v_folds, family = family)
    pred <- as.numeric(predict(fit, newx = as.matrix(W), s = "lambda.min", type = "response"))

  } else if (is.list(method)) {
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
  }

  return(pred)
}
