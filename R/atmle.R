# A-TMLE
atmle <- function(data,
                  S_node,
                  W_node,
                  A_node,
                  Y_node,
                  controls_only,
                  family="gaussian",
                  atmle_pooled = TRUE,
                  var_method = "ic",
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
  theta <- learn_theta(W, A, Y, controls_only, family, nuisance_method)

  if (verbose) print("learning P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A, controls_only, nuisance_method)

  if (verbose) print("learning P(A=1|W)")
  g <- learn_g(S, W, A, g_rct, nuisance_method)

  # learn working model tau for bias
  if (verbose) print("learning E(Y|S,W,A)")
  tau <- learn_tau(S, W, A, Y, Pi, theta, controls_only, working_model)

  # TMLE to target Pi
  if (verbose) print("targeting P(S=1|W,A)")
  for (i in 1:1) {
    # TMLE to target Pi
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi, controls_only)

    # re-learn working model tau with targeted Pi
    tau_star <- learn_tau(S, W, A, Y, Pi_star, theta, controls_only, working_model)
    Pi <- Pi_star
    tau <- tau_star
  }

  psi_pound_est <- NULL
  if (controls_only) {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0)
  } else {
    psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  }

  # estimate pooled-ATE psi_tilde ----------------------------------------------
  psi_tilde <- NULL
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL
  if (atmle_pooled) {
    # use atmle for pooled-ATE
    if (verbose) print("learning E(Y|W)")
    theta_tilde <- learn_theta_tilde(W, Y, family, nuisance_method)

    if (verbose) print("learning psi_tilde")
    psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, "glmnet")

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

  # final estimates ------------------------------------------------------------
  est <- NULL
  lower <- NULL
  upper <- NULL
  psi_pound_lower <- NULL
  psi_pound_upper <- NULL
  psi_tilde_lower <- NULL
  psi_tilde_upper <- NULL

  if (var_method == "ic") {
    # bias parameter
    psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n, controls_only)
    psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
    psi_pound_lower <- psi_pound_est+qnorm(0.025)*psi_pound_se
    psi_pound_upper <- psi_pound_est+qnorm(0.975)*psi_pound_se

    # pooled-ATE parameter
    psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
    psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se

    # RCT-ATE
    est <- psi_tilde_est - psi_pound_est
    eic <- psi_tilde_eic - psi_pound_eic
    se <- sqrt(var(eic, na.rm = TRUE)/n)
    lower <- est+qnorm(0.025)*se
    upper <- est+qnorm(0.975)*se
  } else if (var_method == "bootstrap") {
    # bias parameter
    psi_pound_se <- bootstrap_psi_pound(tau, W, Pi)
    psi_pound_lower <- psi_pound_est+qnorm(0.025)*psi_pound_se
    psi_pound_upper <- psi_pound_est+qnorm(0.975)*psi_pound_se

    # pooled-ATE parameter (use ic-based for now)
    psi_tilde_se <- NULL
    #if (atmle_pooled) {
    #  psi_tilde_se <- bootstrap_psi_tilde(W, psi_tilde)
    #  psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    #  psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se
    #} else {
    psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
    psi_tilde_lower <- psi_tilde_est+qnorm(0.025)*psi_tilde_se
    psi_tilde_upper <- psi_tilde_est+qnorm(0.975)*psi_tilde_se
    #}

    # RCT-ATE
    est <- psi_tilde_est - psi_pound_est
    se <- sqrt(psi_pound_se^2+psi_tilde_se^2)
    lower <- est+qnorm(0.025)*se
    upper <- est+qnorm(0.975)*se
  }

  return(list(est = est,
              lower = lower,
              upper = upper,
              psi_pound_est = psi_pound_est,
              psi_pound_lower = psi_pound_lower,
              psi_pound_upper = psi_pound_upper,
              psi_tilde_est = psi_tilde_est,
              psi_tilde_lower = psi_tilde_lower,
              psi_tilde_upper = psi_tilde_upper))
}
