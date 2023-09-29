# nonparametric estimator, pooling across data
nonparametric <- function(data,
                          S_node,
                          W_node,
                          A_node,
                          Y_node,
                          g_rct,
                          nuisance_method="glm",
                          working_model="lasso",
                          verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)

  # estimate nuisance parts
  if (verbose) print("learning E(Y|S,W,A)")
  Q <- learn_Q_SWA(S, W, A, Y, method = nuisance_method)

  if (verbose) print("learning P(A=1|S,W)")
  g <- learn_g_SW(S, W, A, g_rct, method = nuisance_method)

  if (verbose) print("learning P(S=1|W)")
  Pi <- learn_S_W(S, W, method = nuisance_method)

  # target Q
  Q_star <- target_Q(S, W, A, Y, Pi, g, Q)
  Q_star$pred <- (max(Y)-min(Y))*Q_star$pred+min(Y)
  Q_star$S1A1 <- (max(Y)-min(Y))*Q_star$S1A1+min(Y)
  Q_star$S1A0 <- (max(Y)-min(Y))*Q_star$S1A0+min(Y)

  Q <- Q_star

  # estimates
  psi_est <- mean(Q_star$S1A1-Q_star$S1A0)
  psi_eic <- get_eic_psi_nonparametric(Q_star, Pi, g, S, A, Y, psi_est)
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se
  var_eic <- var(psi_eic)

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper,
              var_eic = var_eic))
}
