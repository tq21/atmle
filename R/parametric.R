# plug-in relaxed HAL psi_pound
# tmle psi_tilde
parametric <- function(data,
                       g_rct=0.67) {
  S <- data$S
  W <- as.matrix(data[, c("W1", "W2", "W3", "W4")])
  A <- data$A
  Y <- data$Y

  n <- nrow(data)

  # estimate bias psi_pound ----------------------------------------------------
  # learn relevant parts
  print("learning E(Y|W)")
  theta_pred <- learn_theta(W, A, Y)
  print("learning P(S=1|W,A)")
  Pi_pred <- learn_Pi(S, W, A)
  print("learning P(A|W)")
  g_pred <- learn_g(S, W, A, g_rct)
  print("learning working model: tau(W,A)=E(Y|S=1,W,A)-E(Y|S=0,W,A)")
  tau_pred <- learn_tau_parametric(S, W, A, Y)

  # plug-in estimates
  psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0-(1-Pi_pred$A1)*tau_pred$A1)
  psi_pound_eic <- get_eic_psi_pound_parametric(Pi_pred, tau_pred, g_pred, psi_pound_est, S, A, Y, n)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  Q_pred <- learn_Q(W, A, Y)
  Q_star <- tmle(Y = Y, A = A, W = W, g1W = g_pred,
                 Q = as.matrix(data.frame(Q_pred$A1, Q_pred$A0)),
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

set.seed(21411)
data <- generate_data(N=1000, p_rct=0.67, bA=0.7)
tictoc::tic()
res <- parametric(data, g_rct = 0.67)
tictoc::toc()
