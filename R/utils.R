
make_counter_design_matrix <- function(basis_list, X_counterfactual, X_unpenalized = NULL) {
  if (!is.null(X_unpenalized)) {
    if (is_empty(basis_list)) {
      # intercept + main terms
      return(cbind(1, X_unpenalized))
    } else {
      # intercept + HAL bases + main terms
      return(cbind(1, cbind(as.matrix(hal9001::make_design_matrix(X_counterfactual, basis_list)), X_unpenalized)))
    }
  } else {
    if (is_empty(basis_list)) {
      # intercept
      return(matrix(rep(1, nrow(X_counterfactual))))
    } else {
      # intercept + HAL bases
      return(cbind(1, as.matrix(hal9001::make_design_matrix(X_counterfactual, basis_list))))
    }
  }
}

to_prob <- function(pred) {
  return(1 / (1 + exp(-pred)))
}

bound <- function(X) {
  X_max <- max(X, na.rm = TRUE)
  X_min <- min(X, na.rm = TRUE)

  return((X-X_min)/(X_max-X_min))
}
