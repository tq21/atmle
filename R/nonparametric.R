# nonparametric, only uses RCT data
nonparametric <- function(data,
                          S_node,
                          W_node,
                          A_node,
                          Y_node,
                          g_rct,
                          verbose=TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[S == 1, W_node]
  A <- data[S == 1, A_node]
  Y <- data[S == 1, Y_node]
  n <- nrow(data)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  # use TMLE to estimate psi_tilde
  Q <- learn_Q(W, A, Y)
  Q_star <- tmle(Y = Y, A = A, W = W, g1W = rep(g_rct, sum(S == 1)),
                 Q = as.matrix(data.frame(Q$A1, Q$A0)),
                 family = "gaussian")
  psi_est <- Q_star$estimates$ATE$psi
  psi_eic <- Q_star$estimates$IC$IC.ATE
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper))
}

# set.seed(1234)
# data <- generate_data(N=500, p_rct=0.67, bA=0.5)
# res <- tmle_nonparametric(data,
#                           S_node = 1,
#                           W_node = c(2, 3, 4, 5),
#                           A_node = 6,
#                           Y_node = 7,
#                           g_rct=0.67)
