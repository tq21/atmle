#' @title Estimate the conditional hazard function for event (failure/censoring)
#'
#' @param data_long A `data.table` of long format data.
#' @param X A `character` vector of covariate names.
#' @param T_tilde A vector of observed failure/censored times.
#' @param event A `character` of event variable name.
#' @param method A `character` of method to use for hazard estimation.
#' @param folds A list of folds for cross-fitting.
#' @param counter_var A `character` of counterfactual variable name, if any.
#'
#' @export
learn_hazard <- function(data_long,
                         X,
                         T_tilde,
                         event,
                         method,
                         folds,
                         counter_var = NULL) {

  if (!is.null(counter_var)) {
    # make counterfactual prediction data
    data_long_A1 <- copy(data_long)
    data_long_A1 <- data_long_A1[, (counter_var) := 1]
    data_long_A0 <- copy(data_long)
    data_long_A0 <- data_long_A0[, (counter_var) := 0]
  }

  # make training data
  data_train <- copy(data_long)
  data_train <- data_train[(T_tilde) >= t]
  cov_names <- c(X, "t")

  if (method == "glm") {
    # fit hazard regression using glm
    fit <- glm(data_train[[event]] ~ .,
               data = data_train[, ..cov_names],
               family = "binomial")
    pred <- predict(fit, newdata = data_long[, ..cov_names], type = "response")

    if (!is.null(counter_var)) {
      A1 <- predict(fit, newdata = data_long_A1[, ..cov_names], type = "response")
      A0 <- predict(fit, newdata = data_long_A0[, ..cov_names], type = "response")
    }

  } else if (method == "HAL") {
    # fit hazard regression using HAL
    fit <- fit_hal(X = data_train[, ..cov_names],
                   Y = data_train[[event]],
                   id = data_train[["id"]],
                   smoothness_orders = 1,
                   max_degree = 3,
                   family = "binomial",
                   return_x_basis = FALSE)

    # predict hazard
    pred <- predict(fit, new_data = data_long[, ..cov_names], type = "response")

    if (!is.null(counter_var)) {
      A1 <- predict(fit, new_data = data_long_A1[, ..cov_names], type = "response")
      A0 <- predict(fit, new_data = data_long_A0[, ..cov_names], type = "response")
    }
  }

  if (is.null(counter_var)) {
    return(as.numeric(pred))
  } else {
    return(list(pred = as.numeric(pred),
                A1 = as.numeric(A1),
                A0 = as.numeric(A0)))
  }
}
