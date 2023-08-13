atmle <- function(data, g_rct, max_iter=1, eps=1e-5) {
  S <- data$S
  W <- as.matrix(data[, c("W1", "W2", "W3", "W4")])
  A <- data$A
  Y <- data$Y

  n <- nrow(data)

  # estimate bias psi_pound ----------------------------------------------------
  # learn relevant parts
  print("learning E(Y|W,A)")
  theta_pred <- learn_theta(W, A, Y)
  print("learning P(S=1|W,A)")
  Pi_pred <- learn_Pi(S, W, A)
  print("learning P(A|W)")
  g_pred <- learn_g(S, W, A, g_rct)
  print("learning E(Y|S,W,A)")
  tau_pred <- learn_tau(S, W, A, Y, Pi_pred$pred, theta_pred)

  # TMLE target Pi
  cur_iter <- 1
  tol <- Inf
  while (cur_iter <= max_iter && tol > eps) {
    print(cur_iter)

    Pi_star <- Pi_tmle(S, W, A, g_pred, tau_pred, Pi_pred)
    tau_star <- learn_tau(S, W, A, Y, Pi_star$pred, theta_pred)
    cur_eic <- get_eic_Pi(g_pred, tau_star, Pi_star, S, A)
    tol <- abs(mean(cur_eic) - 0)
    cur_iter <- cur_iter + 1

    print(tol)
  }

  psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0-(1-Pi_star$A1)*tau_star$A1)
  psi_pound_eic <- get_eic_psi_pound(Pi_star, tau_star, g_pred, theta_pred, psi_pound_est, S, A, Y, n)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  print("learning E(Y|W)")
  theta_tilde_pred <- learn_theta_tilde(W, Y)
  print("learning psi_tilde")
  psi_tilde <- learn_psi_tilde(W, A, Y, g_pred, theta_tilde_pred)
  psi_tilde_est <- mean(psi_tilde$pred)
  psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g_pred, theta_pred, Y, A)

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

data <- generate_data(N=500, p_rct=0.67, bA=1.5)
res <- atmle(data, g_rct=0.67, max_iter=1, eps=1e-5)
