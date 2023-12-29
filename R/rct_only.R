# estimator that only uses RCT data
rct_only <- function(data,
                     S_node,
                     W_node,
                     A_node,
                     Y_node,
                     g_rct,
                     family,
                     nuisance_method="glm",
                     verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[S == 1, W_node]
  A <- data[S == 1, A_node]
  Y <- data[S == 1, Y_node]
  n <- nrow(W)

  # estimate ATE using TMLE
  Q <- learn_Q(W, A, Y, nuisance_method, 5, family, c(-Inf, Inf))
  Q_star <- tmle(Y = Y, A = A, W = W, g1W = rep(g_rct, sum(S == 1)),
                 Q = as.matrix(data.frame(Q$A1, Q$A0)),
                 family = family)
  psi_est <- Q_star$estimates$ATE$psi
  psi_eic <- Q_star$estimates$IC$IC.ATE
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper))
}
