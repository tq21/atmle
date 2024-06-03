learn_Q_bar_r_working <- function(W,
                                  Q_bar_r,
                                  G_bar) {

  # TODO: here, we might be able to add bound by using sigmoid link
  fit <- fit_hal(X = W,
                 Y = Q_bar_r,
                 smoothness_orders = 0,
                 max_degree = 3,
                 family = "gaussian",
                 weights = 1/G_bar$integrate_A,
                 return_x_basis = TRUE)

  coefs <- fit$coefs[which(fit$coefs != 0)]
  phi <- as.matrix(cbind(1, fit$x_basis[, which(fit$coefs[-1] != 0), drop = FALSE]))
  pred <- as.numeric(phi %*% coefs)

  return(list(coefs = coefs,
              phi = phi,
              pred = pred))
}
