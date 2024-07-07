# nonparametric estimator, pooling across data
nonparametric <- function(data,
                          S,
                          W,
                          A,
                          Y,
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
                          var_method = "ic",
                          max_iter = 1,
                          v_folds = 5,
                          verbose = TRUE) {

  # define nodes ---------------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  delta <- as.numeric(!is.na(Y)) # missingness indicator
  n <- nrow(data) # sample size

  # estimate nuisance parts
  if (verbose) print("learning E(Y|S,W,A)")
  Q <- learn_Q_S1(
    S = S,
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    method = Q_method
  )

  if (verbose) print("learning P(A=1|S,W)")
  g <- learn_g(W = W[S == 1, ],
               A = A[S == 1],
               method = g_method,
               folds = 1:5,
               g_bounds = c(0, 1),
               cross_fit_nuisance = FALSE)

  if (verbose) print("learning P(S=1|W)")
  Pi <- learn_S_W(S, W, Pi_method)

  # estimate missing mechanism
  # g_delta <- learn_g_delta(
  #   W = W,
  #   A = A,
  #   delta = delta,
  #   method = g_method,
  #   folds = c(1, 2, 3, 4, 5),
  #   g_bounds = c(0, 1)
  # )

  # censoring weights
  g_delta <- list(pred = rep(1, n))
  weights <- delta / g_delta$pred

  # target Q
  Q_star <- target_Q(S, W, A, Y, Pi, g, Q, delta, g_delta)
  Q <- Q_star

  # estimates
  psi_est <- mean(Q$S1A1 - Q$S1A0)
  psi_eic <- get_eic_psi_nonparametric(Q, Pi, g, S, A, Y, psi_est, weights)
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE) / n)
  psi_ci_lower <- psi_est - 1.96 * psi_se
  psi_ci_upper <- psi_est + 1.96 * psi_se

  return(list(
    est = psi_est,
    lower = psi_ci_lower,
    upper = psi_ci_upper
  ))
}
