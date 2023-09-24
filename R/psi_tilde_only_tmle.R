# TMLE for psi_tilde
psi_tilde_only_tmle <- function(data,
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

  # estimate pooled ATE psi_tilde ----------------------------------------------
  if (verbose) print("learning P(A|W)")
  g <- learn_g(S, W, A, p_rct, nuisance_method)

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

  return(list(est = psi_tilde_est,
              lower = psi_tilde_ci_lower,
              upper = psi_tilde_ci_upper))
}
