
fit_relaxed_hal <- function(X, Y,
                            family,
                            weights = NULL,
                            X_no_basis = NULL,
                            X_unpenalized = NULL) {
  x_basis <- NULL
  basis_list <- NULL

  if (!is.null(X_no_basis)) {
    # include main terms in HAL
    hal_fit <- fit_hal(X = X,
                       Y = Y,
                       X_no_basis = X_no_basis,
                       family = family,
                       weights = weights,
                       max_degree = 3,
                       smoothness_orders = 0)
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
                       #num_knots = c(100, 50, 50))
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
