learn_Q_bar_r_working <- function(W,
                                  Q_bar_r,
                                  G_bar,
                                  weights) {

  # TODO: here, we might be able to add bound by using sigmoid link
  fit <- fit_hal(X = W,
                 Y = Q_bar_r,
                 smoothness_orders = 0,
                 max_degree = 3,
                 family = "gaussian",
                 weights = weights,
                 return_x_basis = TRUE)

  # coefs <- fit$coefs[which(fit$coefs != 0)]
  # pred <- as.numeric(phi %*% coefs)

  phi <- as.matrix(fit$x_basis[, which(fit$coefs[-1] != 0), drop = FALSE])

  # fit unpenalized
  relaxed_fit <- glm(Q_bar_r ~ phi, family = "gaussian", weights = weights)
  coefs <- as.numeric(coef(relaxed_fit))
  coefs[is.na(coefs)] <- 0
  phi <- as.matrix(cbind(1, phi[, which(coefs[-1] != 0), drop = FALSE]))
  new_coefs <- coefs[which(coefs != 0)]
  pred <- as.numeric(phi %*% matrix(new_coefs))

  return(list(coefs = new_coefs,
              phi = phi,
              pred = pred))
}
