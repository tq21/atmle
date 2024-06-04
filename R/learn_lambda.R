#' @title Estimate the conditional failure-time hazard function
#'
#' @param W A matrix of covariates.
#' @param A A vector of treatment assignments.
#' @param T_tilde A vector of observed failure times.
#' @param folds A list of folds for cross-fitting.
#'
#' @export
learn_lambda <- function(W,
                         A,
                         T_tilde,
                         Delta,
                         method,
                         folds) {

  # prepare repeated data for pooled logistic regression
  rep_data_both <- make_rep_data(W = cbind(W, A = A), T_tilde = T_tilde, Delta = Delta)
  cov_names <- c(names(W), "A")

  # prepare data for predictions
  data_pred <- rep_data_both$rep_data
  data_A1 <- copy(data_pred); data_A1 <- data_A1[, A := 1]
  data_A0 <- copy(data_pred); data_A0 <- data_A0[, A := 0]

  if (method == "HAL") {

    # fit hazard regression
    fit <- fit_hal(X = rep_data_both$rep_data_small[, ..cov_names],
                   Y = rep_data_both$rep_data_small$T_tilde_t,
                   id = rep_data_both$rep_data_small$id,
                   smoothness_orders = 0,
                   family = "binomial",
                   return_x_basis = TRUE)

    # hazard predictions
    A1 <- predict(fit, new_data = data_A1, type = "response")
    A0 <- predict(fit, new_data = data_A0, type = "response")
  } else if (method == "glm") {

    # fit logistic regression
    var_names <- c(cov_names, "T_tilde_t")
    fit <- glm(T_tilde_t ~ .,
               data = rep_data_both$rep_data_small[, ..var_names],
               family = "binomial")

    # hazard predictions
    A1 <- as.numeric(predict(fit, newdata = data_A1, type = "response"))
    A0 <- as.numeric(predict(fit, newdata = data_A0, type = "response"))
  }

  data_pred <- data_pred[, `:=` (lambda_A1 = A1, lambda_A0 = A0)]

  return(data_pred)
}
