#' @title Learn nuisance function: trial enrollment mechanism
#'
#' @description Function to learn the trial enrollment mechanism
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
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
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
#' @param method Learning method. \code{"glm"} for main-term linear model,
#' \code{"glmnet"} for lasso, or a \code{list} of \code{sl3} learners for
#' super learner-based estimation.
#' @param v_folds A numeric of number of folds for cross-validation
#' (when necessary).
#' @param Pi_bounds A numeric vector of lower and upper bounds for the
#' trial enrollment probabilities. The first element is the lower bound,
#' and the second element is the upper bound.
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{The estimated trial enrollment probabilities;}
#' \item{A1}{The estimated trial enrollment probabilities under treatment;}
#' \item{A0}{The estimated trial enrollment probabilities under control.}
learn_Pi <- function(S,
                     W,
                     A,
                     controls_only,
                     method,
                     v_folds,
                     Pi_bounds) {

  if (is.character(method) && method == "sl3") {
    method <- get_default_sl3_learners("binomial")
  }

  pred <- numeric(length(A))
  A1 <- rep(1, length(A))
  A0 <- numeric(length(A))

  if (is.list(method)) {
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
    A0 <- .bound(fit_Pi$predict(task_Pi_A0), Pi_bounds)
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1 # no treated in external
    } else {
      task_Pi_A1 <- sl3_Task$new(data = data.table(W, A = 1, S = S),
                                 covariates = c(colnames(W), "A"),
                                 outcome = "S", outcome_type = "binomial")
      A1 <- .bound(fit_Pi$predict(task_Pi_A1), Pi_bounds)
      pred[A == 1] <- A1[A == 1]
    }
  } else if (method == "glm") {
    # A = 0
    fit_A0 <- glm(S[A == 0] ~., data = data.frame(W[A == 0,]),
                  family = "binomial")
    A0 <- .bound(as.numeric(predict(fit_A0, newdata = data.frame(W),
                                    type = "response")), Pi_bounds)
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1 # no treated in external
    } else {
      # A = 1
      fit_A1 <- glm(S[A == 1] ~., data = data.frame(W[A == 1,]),
                    family = "binomial")
      A1 <- .bound(as.numeric(predict(fit_A1, newdata = data.frame(W),
                                      type = "response")), Pi_bounds)
      pred[A == 1] <- A1[A == 1]
    }

  } else if (method == "glmnet") {
    # A = 0
    fit_A0 <- cv.glmnet(x = as.matrix(W[A == 0, ]), y = S[A == 0],
                        keep = TRUE, alpha = 1, nfolds = v_folds,
                        family = "binomial")
    A0 <- .bound(as.numeric(predict(fit_A0, newx = as.matrix(W),
                                    s = "lambda.min",
                                    type = "response")), Pi_bounds)
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1
    } else {
      # A = 1
      fit_A1 <- cv.glmnet(x = as.matrix(W[A == 1, ]), y = S[A == 1],
                          keep = TRUE, alpha = 1, nfolds = v_folds,
                          family = "binomial")
      A1 <- .bound(as.numeric(predict(fit_A1, newx = as.matrix(W),
                                      s = "lambda.min",
                                      type = "response")), Pi_bounds)
      pred[A == 1] <- A1[A == 1]
    }

  } else if (method == "empirical") {
    # estimate Pi using its empirical distribution, not recommended
    # for testing purposes only, only use if S independent of W
    # A == 0
    A0 <- rep(mean(S[A == 0]), length(A))
    pred[A == 0] <- A0[A == 0]

    if (controls_only) {
      pred[A == 1] <- 1
    } else {
      A1 <- rep(mean(S[A == 1]), length(A))
      pred[A == 1] <- A1[A == 1]
    }
  } else {
    stop("Invalid method. Must be one of 'glm', 'glmnet', or 'sl3', or a
         list of sl3 learners.")
  }

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}
