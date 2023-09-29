fit_relaxed_hal <- function(X, Y, family, weights = NULL, add_main_terms = FALSE) {
  # HAL fit
  x_basis <- NULL
  basis_list <- NULL
  if (add_main_terms) {
    # include main terms in HAL
    hal_fit <- fit_hal(X = X,
                       Y = Y,
                       family = family,
                       weights = weights,
                       max_degree = 3,
                       smoothness_orders = 0,
                       X_unpenalized = X)
    basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
    basis_list <- basis_list[1:(length(basis_list)-ncol(X))] # remove main terms

    if (length(basis_list) == 1 & is.null(basis_list[[1]])) {
      # basis list empty, intercept + main terms
      x_basis <- as.matrix(cbind(1, X))
    } else {
      # basis list non-empty, intercept + HAL bases + main terms
      x_basis <- cbind(1, cbind(as.matrix(make_design_matrix(X, basis_list)), X))
    }
  } else {
    hal_fit <- fit_hal(X = X,
                       Y = Y,
                       family = family,
                       weights = weights,
                       max_degree = 3,
                       smoothness_orders = 0)
    basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]

    # basis list empty
    if (length(basis_list) == 0) {
      x_basis <- matrix(rep(1, length(Y)))
    } else {
      x_basis <- cbind(1, as.matrix(make_design_matrix(X, basis_list)))
    }
  }

  # relaxed HAL fit
  beta <- NULL
  pred <- NULL
  if (family == "binomial") {
    hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = binomial(), weights = weights, intercept = FALSE)
    beta <- coef(hal_relaxed_fit)
    beta[is.na(beta)] <- 0
    pred <- x_basis %*% beta
    pred <- as.vector(1 / (1 + exp(-pred)))
  } else if (family == "gaussian") {
    hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = gaussian(), weights = weights, intercept = FALSE)
    beta <- coef(hal_relaxed_fit)
    beta[is.na(beta)] <- 0
    pred <- as.vector(x_basis %*% beta)
  }

  return(list(beta = beta, basis_list = basis_list, pred = pred))
}

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
