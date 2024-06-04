#' @title Make repeated measures data for pooled logistic regression
#' for estimating conditional fT_tilde_ilure-time hazard function.
#'
#' @param W A matrix of covariates.
#' @param T_tilde A vector of time-points.
#'
#' @export
make_rep_data <- function(W, T_tilde, Delta) {

  tau <- length(unique(T_tilde))

  # make repeated data for pooled logistic regression
  rep_data <- data.table(W)
  rep_data <- rep_data[rep(1:.N, each = tau)]
  rep_data <- rep_data[, `:=` (id = rep(seq(length(T_tilde)), each = tau),
                               t = rep(seq(tau), n),
                               T_tilde = rep(T_tilde, each = tau),
                               Delta = rep(Delta, each = tau))]
  rep_data <- rep_data[, T_tilde_t := as.numeric(T_tilde == t & Delta == 1)]
  rep_data_small <- copy(rep_data)
  rep_data_small <- rep_data_small[, filter_rows(.SD, "T_tilde_t"), by = id]

  return(list(rep_data = rep_data,
              rep_data_small = rep_data_small))
}
