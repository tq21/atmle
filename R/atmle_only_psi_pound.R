# A-TMLE only psi_pound, TMLE psi_tilde
atmle_tmle <- function(data,
                       S_node,
                       W_node,
                       A_node,
                       Y_node,
                       g_rct,
                       target_Pi,
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
  theta <- learn_theta(W, A, Y)

  if (verbose) print("learning P(S=1|W,A)")
  Pi <- learn_Pi(S, W, A)

  if (verbose) print("learning P(A|W)")
  g <- learn_g(S, W, A, g_rct)

  # learn initial estimate of working model tau
  if (verbose) print("learning E(Y|S,W,A)")
  tau <- learn_tau(S, W, A, Y, Pi$pred, theta)

  psi_pound_est <- NULL
  psi_pound_eic <- NULL
  if (target_Pi) {
    # TMLE to target Pi
    if (verbose) print("targeting P(S=1|W,A)")
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi)

    # re-learn working model tau with targeted Pi
    tau_star <- learn_tau(S, W, A, Y, Pi_star$pred, theta)

    psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0-(1-Pi_star$A1)*tau_star$A1)
    psi_pound_eic <- get_eic_psi_pound(Pi_star, tau_star, g, theta, psi_pound_est, S, A, Y, n)
  } else {
    # No targeting of Pi, for comparison only
    psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
    psi_pound_eic <- get_eic_psi_pound(Pi, tau, g, theta, psi_pound_est, S, A, Y, n)
  }

  # estimate pooled ATE psi_tilde ----------------------------------------------
  # use TMLE to estimate psi_tilde
  Q <- learn_Q(W, A, Y)
  Q_star <- tmle(Y = Y, A = A, W = W, g1W = g,
                 Q = as.matrix(data.frame(Q$A1, Q$A0)),
                 family = "gaussian")
  psi_tilde_est <- Q_star$estimates$ATE$psi
  psi_tilde_eic <- Q_star$estimates$IC$IC.ATE

  # estimate psi ---------------------------------------------------------------
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_pound_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_pound_se
  psi_ci_upper <- psi_est+1.96*psi_pound_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper))
}

# data <- generate_data(N=500, p_rct=0.67, bA=0.5)
# res <- atmle_tmle(data,
#                   S_node = 1,
#                   W_node = c(2, 3, 4, 5),
#                   A_node = 6,
#                   Y_node = 7,
#                   target_Pi = TRUE,
#                   g_rct=0.67)
