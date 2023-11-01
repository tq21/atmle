# nonparametric estimator, pooling across data
nonparametric <- function(data,
                          S_node,
                          W_node,
                          A_node,
                          Y_node,
                          controls_only,
                          family = "gaussian",
                          atmle_pooled = TRUE,
                          theta_method = "glm",
                          Pi_method = "glm",
                          g_method = "glm",
                          theta_tilde_method = "glm",
                          Q_method = "glm",
                          bias_working_model = "glmnet",
                          pooled_working_model = "glmnet",
                          g_rct,
                          var_method = "ic",
                          max_iter = 1,
                          v_folds = 5,
                          verbose = TRUE) {

  # define nodes ---------------------------------------------------------------
  S <- data[, S_node] # study indicator
  W <- data[, W_node] # covariates
  A <- data[, A_node] # treatment
  Y <- data[, Y_node] # outcome
  n <- nrow(data) # sample size

  # estimate nuisance parts
  if (verbose) print("learning E(Y|S,W,A)")
  Q <- learn_Q_S1(S, W, A, Y, Q_method)

  if (verbose) print("learning P(A=1|S,W)")
  g <- rep(g_rct, n)

  if (verbose) print("learning P(S=1|W)")
  Pi <- learn_S_W(S, W, Pi_method)

  # target Q
  Q_star <- target_Q(S, W, A, Y, Pi, g, Q)
  Q <- Q_star

  # estimates
  psi_est <- mean(Q$S1A1-Q$S1A0)
  psi_eic <- get_eic_psi_nonparametric(Q, Pi, g, S, A, Y, psi_est)
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(est = psi_est,
              lower = psi_ci_lower,
              upper = psi_ci_upper))
}
