
make_counter_design_matrix <- function(basis_list, X_counterfactual, add_main_terms = FALSE) {
  if (add_main_terms) {
    if (length(basis_list) == 1 & is.null(basis_list[[1]])) {
      # intercept + main terms
      return(cbind(1, X_counterfactual))
    } else {
      # intercept + HAL bases + main terms
      return(cbind(1, cbind(as.matrix(hal9001::make_design_matrix(X_counterfactual, basis_list)), X_counterfactual)))
    }
  } else {
    if (length(basis_list) == 0) {
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
