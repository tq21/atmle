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
                         method,
                         folds) {

  if (method == "HAL") {

    # prepare repeated data for pooled logistic regression
    rep_data <- make_rep_data(W = cbind(W, A = A), T_tilde = T_tilde)

    # fit hazard regression
    fit <- fit_hal(X = rep_data$data,
                   Y = rep_data$ind,
                   id = rep_data$id,
                   smoothness_orders = 0,
                   family = "binomial",
                   return_x_basis = TRUE)

    # make data for predictions
    unique_T_tilde <- sort(unique(T_tilde))
    data_pred <- do.call(rbind, replicate(length(unique_T_tilde), cbind(W, A = A), simplify = FALSE))
    data_pred <- cbind(T_tilde_i = rep(unique_T_tilde, each = length(A)), data_pred)
    data_A1 <- data_pred; data_A1$A <- 1
    data_A0 <- data_pred; data_A0$A <- 0

    # hazard predictions
    A1 <- predict(fit, new_data = data_A1, type = "response")
    A0 <- predict(fit, new_data = data_A0, type = "response")
    id <- rep(seq(length(A)), length(unique_T_tilde))
    T_tilde_i <- data_pred$T_tilde_i
  }

  return(list(A1 = A1,
              A0 = A0,
              id = id,
              T_tilde_i = T_tilde_i))
}
