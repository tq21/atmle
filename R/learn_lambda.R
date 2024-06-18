#' @title Estimate the conditional failure-time hazard function
#'
#' @param W A matrix of covariates.
#' @param A A vector of treatment assignments.
#' @param T_tilde A vector of observed failure times.
#' @param folds A list of folds for cross-fitting.
#'
#' @export
learn_lambda <- function(data_long,
                         W,
                         A,
                         T_tilde,
                         method,
                         folds) {

  # make counterfactual prediction data
  data_long_A1 <- copy(data_long)
  data_long_A1 <- data_long_A1[, A := 1]
  data_long_A0 <- copy(data_long)
  data_long_A0 <- data_long_A0[, A := 0]

  # make training data
  data_train <- copy(data_long)
  data_train <- data_train[(T_tilde) >= t]
  cov_names <- c(W, A, "t")

  if (method == "glm") {
    # fit hazard regression using glm
    fit <- glm(data_train[["T_t"]] ~ .,
               data = data_train[, ..cov_names],
               family = "binomial")
    pred <- predict(fit, newdata = data_long[, ..cov_names], type = "response")
    A1 <- predict(fit, newdata = data_long_A1[, ..cov_names], type = "response")
    A0 <- predict(fit, newdata = data_long_A0[, ..cov_names], type = "response")

  } else if (method == "HAL") {
    # fit hazard regression using HAL
    fit <- fit_hal(X = data_train[, ..cov_names],
                   Y = data_train[["T_t"]],
                   id = data_train[["id"]],
                   smoothness_orders = 1,
                   max_degree = 3,
                   family = "binomial",
                   return_x_basis = FALSE)

    # predict hazard
    pred <- predict(fit, new_data = data_long[, ..cov_names], type = "response")
    A1 <- predict(fit, new_data = data_long_A1[, ..cov_names], type = "response")
    A0 <- predict(fit, new_data = data_long_A0[, ..cov_names], type = "response")
  }

  return(list(pred = as.numeric(pred),
              A1 = as.numeric(A1),
              A0 = as.numeric(A0)))
}
