# A-TMLE both psi_pound and psi_tilde
atmle <- function(data,
                  S_node,
                  W_node,
                  A_node,
                  Y_node,
                  nuisance_method="lasso",
                  working_model="lasso",
                  p_rct=0.5,
                  verbose=TRUE,
                  transform=TRUE) {

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
  tau <- NULL
  if (transform) {
    tau <- learn_tau(S, W, A, Y, Pi, theta, working_model)
  } else {
    tau <- learn_tau_test(S, W, A, Y, Pi, theta, working_model)
  }
  psi_pound_est <- NULL
  psi_pound_eic <- NULL

  # TMLE to target Pi
  #if (verbose) print("targeting P(S=1|W,A)")
  #for (i in 1:1) {
  #  # TMLE to target Pi
  #  Pi_star <- Pi_tmle(S, W, A, g, tau, Pi)
#
  #  # re-learn working model tau with targeted Pi
  #  tau_star <- NULL
  #  if (transform) {
  #    tau_star <- learn_tau(S, W, A, Y, Pi_star, theta, working_model)
  #  } else {
  #    tau_star <- learn_tau_test(S, W, A, Y, Pi_star, theta, working_model)
  #  }
#
  #  Pi <- Pi_star
  #  tau <- tau_star
  #}

  psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n)
  print("psi pound eic: " %+% mean(psi_pound_eic))

  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE))
  psi_pound_ci_lower <- psi_pound_est-1.96*psi_pound_se/sqrt(n)
  psi_pound_ci_upper <- psi_pound_est+1.96*psi_pound_se/sqrt(n)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning E(Y|W)")
  theta_tilde <- learn_theta_tilde(W, Y, nuisance_method)

  # learn psi_tilde using R-loss
  if (verbose) print("learning psi_tilde")
  psi_tilde <- NULL
  if (transform) {
    psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, "lasso")
  } else {
    psi_tilde <- learn_psi_tilde_test(W, A, Y, g, theta_tilde, "lasso")
  }
  psi_tilde_est <- mean(psi_tilde$pred)
  psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g, theta_tilde, Y, A, n)

  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE))
  psi_tilde_ci_lower <- psi_tilde_est-1.96*psi_tilde_se/sqrt(n)
  psi_tilde_ci_upper <- psi_tilde_est+1.96*psi_tilde_se/sqrt(n)
  print("psi tilde eic: " %+% mean(psi_tilde_eic))

  # estimate psi ---------------------------------------------------------------
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se
  var_eic <- var(psi_eic)

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper,
              psi_pound_est = psi_pound_est,
              psi_pound_lower = psi_pound_ci_lower,
              psi_pound_upper = psi_pound_ci_upper,
              psi_tilde_est = psi_tilde_est,
              psi_tilde_lower = psi_tilde_ci_lower,
              psi_tilde_upper = psi_tilde_ci_upper,
              var_eic = var_eic))
}
