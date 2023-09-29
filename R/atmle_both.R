# A-TMLE
atmle <- function(data,
                  S_node,
                  W_node,
                  A_node,
                  Y_node,
                  atmle_pooled = TRUE,
                  nuisance_method="glmnet",
                  working_model="glmnet",
                  g_rct=0.5,
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
  theta <- learn_theta(W, A, Y, family = "gaussian", nuisance_method)

  if (verbose) print("learning P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A, nuisance_method)

  if (verbose) print("learning P(A|W)")
  g <- learn_g(S, W, A, g_rct, nuisance_method)

  # learn working model tau for bias
  if (verbose) print("learning E(Y|S,W,A)")
  tau <- learn_tau(S, W, A, Y, Pi, theta, working_model)

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
  print("psi pound eic: " %+% mean(psi_pound_eic))

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) print("learning E(Y|W)")
    theta_tilde <- learn_theta_tilde(W, Y, family = "gaussian", nuisance_method)

    if (verbose) print("learning psi_tilde")
    psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, "lasso")

    # estimates
    psi_tilde_est <- mean(psi_tilde$pred)
    psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g, theta_tilde, Y, A, n)
  } else {
    # use regular TMLE for pooled-ATE
    Q <- learn_Q(W, A, Y, method = nuisance_method)
    Q_star <- tmle(Y = Y, A = A, W = W, g1W = g,
                   Q = as.matrix(data.frame(Q$A1, Q$A0)),
                   family = "gaussian")

    # estimates
    psi_tilde_est <- Q_star$estimates$ATE$psi
    psi_tilde_eic <- Q_star$estimates$IC$IC.ATE
  }

  print("psi tilde eic: " %+% mean(psi_tilde_eic))

  # final estimates ------------------------------------------------------------
  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
  psi_pound_ci_lower <- psi_pound_est+qnorm(0.025)*psi_pound_se
  psi_pound_ci_upper <- psi_pound_est+qnorm(0.975)*psi_pound_se
  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
  psi_tilde_ci_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
  psi_tilde_ci_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est+qnorm(0.025)*psi_se
  psi_ci_upper <- psi_est+qnorm(0.975)*psi_se

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
