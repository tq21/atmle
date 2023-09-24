# A-TMLE only psi_pound, TMLE psi_tilde
atmle_tmle <- function(data,
                       S_node,
                       W_node,
                       A_node,
                       Y_node,
                       nuisance_method="lasso",
                       working_model="HAL",
                       p_rct=0.5,
                       verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)

  # estimate bias psi_pound ----------------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning E(Y|W,A)")
  theta <- learn_theta(W, A, Y, nuisance_method)

  if (verbose) print("learning P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A, nuisance_method)

  if (verbose) print("learning P(A|W)")
  g <- learn_g(S, W, A, p_rct, nuisance_method)

  # learn initial estimate of working model tau
  if (verbose) print("learning E(Y|S,W,A)")
  tau <- learn_tau(S, W, A, Y, Pi, theta, working_model)

  psi_pound_est <- NULL
  psi_pound_eic <- NULL

  # TMLE to target Pi
  if (verbose) print("targeting P(S=1|W,A)")
  for (i in 1:1) {
    # TMLE to target Pi
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi)

    # re-learn working model tau with targeted Pi
    tau_star <- learn_tau(S, W, A, Y, Pi_star, theta, working_model)

    Pi <- Pi_star
    tau <- tau_star
  }

  psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n)

  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE))
  psi_pound_ci_lower <- psi_pound_est-1.96*psi_pound_se/sqrt(n)
  psi_pound_ci_upper <- psi_pound_est+1.96*psi_pound_se/sqrt(n)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  # use TMLE to estimate psi_tilde
  Q <- learn_Q(W, A, Y, method = nuisance_method)
  Q_star <- tmle(Y = Y, A = A, W = W, g1W = g,
                 Q = as.matrix(data.frame(Q$A1, Q$A0)),
                 family = "gaussian")
  psi_tilde_est <- Q_star$estimates$ATE$psi
  psi_tilde_eic <- Q_star$estimates$IC$IC.ATE

  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
  psi_tilde_ci_lower <- psi_tilde_est-1.96*psi_tilde_se
  psi_tilde_ci_upper <- psi_tilde_est+1.96*psi_tilde_se

  # estimate psi ---------------------------------------------------------------
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper,
              psi_pound_est = psi_pound_est,
              psi_pound_lower = psi_pound_ci_lower,
              psi_pound_upper = psi_pound_ci_upper,
              psi_tilde_est = psi_tilde_est,
              psi_tilde_lower = psi_tilde_ci_lower,
              psi_tilde_upper = psi_tilde_ci_upper))
}
