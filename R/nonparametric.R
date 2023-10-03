# nonparametric estimator, pooling across data
nonparametric <- function(data,
                          S_node,
                          W_node,
                          A_node,
                          Y_node,
                          g_rct,
                          nuisance_method="glm",
                          working_model="glmnet",
                          verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)

  # estimate nuisance parts
  if (verbose) print("learning E(Y|S,W,A)")
  Q <- learn_Q_S1(S, W, A, Y, method = nuisance_method)

  if (verbose) print("learning P(A=1|S,W)")
  g <- rep(g_rct, n)

  if (verbose) print("learning P(S=1|W)")
  Pi <- learn_S_W(S, W, method = nuisance_method)

  # target Q
  Q_star <- target_Q(S, W, A, Y, Pi, g, Q)
  Q <- Q_star

  # estimates
  psi_est <- mean(Q$S1A1-Q$S1A0)
  psi_eic <- get_eic_psi_nonparametric(Q, Pi, g, S, A, Y, psi_est)
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper))
}
