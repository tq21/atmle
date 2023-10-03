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
learn_g <- function(S, W, A, g_rct, method = "glmnet", v_folds = 5) {
  pred <- NULL

  if (method == "glmnet") {
    fit <- cv.glmnet(x = as.matrix(W), y = A,
                     keep = TRUE, alpha = 1, nfolds = v_folds, family = "binomial")
    pred <- as.numeric(predict(fit, newx = as.matrix(W), s = "lambda.min", type = "response"))
  } else if (method == "glm") {
    fit <- glm(A ~., data = W, family = "binomial")
    pred <- as.numeric(predict(fit, newdata = W, type = "response"))
  } else if (method == "sl3") {
    lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 3, family = "binomial"),
                                 Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0)))
    lrnr_g <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new(loss_loglik_binomial))
    task_g <- sl3_Task$new(data = data.table(W, A = A),
                           covariates = colnames(W),
                           outcome = "A", outcome_type = "binomial")
    fit_g <- lrnr_g$train(task_g)
    pred <- fit_g$predict(task_g)
  }

  return(pred)
}

# learn_g <- function(S, W, A, g_rct, method = "glmnet", v_folds = 5) {
#   pred <- NULL
#
#   if (method == "glmnet") {
#     fit <- cv.glmnet(x = as.matrix(W), y = A,
#                      keep = TRUE, alpha = 1, nfolds = v_folds, family = "binomial")
#     pred <- as.numeric(predict(fit, newx = as.matrix(W), s = "lambda.min", type = "response"))
#   } else if (method == "glm") {
#     # P(S=1|W)
#     fit_s_w <- glm(S ~., data = W, family = "binomial")
#     pred_s_w <- as.numeric(predict(fit_s_w, newdata = W, type = "response"))
#     fit_a_ws0 <- glm(A[S == 0] ~., data = W[S == 0,], family = "binomial")
#     pred_a_ws0 <- as.numeric(predict(fit_a_ws0, newdata = W, type = "response"))
#     pred <- g_rct*pred_s_w+pred_a_ws0*(1-pred_s_w)
#
#     #fit <- glm(A ~., data = W, family = "binomial")
#     #pred <- as.numeric(predict(fit, newdata = W, type = "response"))
#   } else if (method == "sl3") {
#     lrnr_stack <- Stack$new(list(Lrnr_earth$new(degree = 3, family = "binomial"),
#                                  Lrnr_xgboost$new(max_depth = 4, nrounds = 20, verbose = 0)))
#     lrnr_g <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new(loss_loglik_binomial))
#     task_g <- sl3_Task$new(data = data.table(W, A = A),
#                            covariates = colnames(W),
#                            outcome = "A", outcome_type = "binomial")
#     fit_g <- lrnr_g$train(task_g)
#     pred <- fit_g$predict(task_g)
#   }
#
#   return(pred)
# }







