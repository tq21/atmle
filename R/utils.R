
fit_relaxed_hal <- function(X, Y, family, weights=NULL, smoothness_orders=0, num_knots=c(200,100,50)) {

  # fit hal
  hal_fit <- fit_hal(X = X, Y = Y, family = family, max_degree = 3, weights = weights,
                     smoothness_orders = smoothness_orders, num_knots = num_knots)
  basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  x_basis <- cbind(1, as.matrix(hal9001::make_design_matrix(X, basis_list)))

  # relaxed fit
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

make_counter_design_matrix <- function(basis_list, X_counterfactual) {
  return(cbind(1, as.matrix(hal9001::make_design_matrix(X_counterfactual, basis_list))))
}

to_prob <- function(pred) {
  return(1 / (1 + exp(-pred)))
}
