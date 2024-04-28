procova <- function(data,
                    S_node,
                    W_node,
                    A_node,
                    Y_node,
                    controls_only,
                    family,
                    g_rct,
                    Q_method = "glmnet",
                    g_method = "glmnet",
                    g_delta_method = "glmnet",
                    v_folds = 5,
                    g_bounds = c(0.01, 0.99),
                    target_gwt = TRUE,
                    verbose = TRUE) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # define nodes ---------------------------------------------------------------
  S <- data[, S_node] # study indicator
  W <- data[, W_node] # covariates
  A <- data[, A_node] # treatment
  Y <- data[, Y_node] # outcome
  delta <- as.integer(!is.na(Y)) # missingness indicator
  n <- nrow(data) # sample size

  # estimate additional score from external data
  procova_score <- learn_procova_score(
    S = S,
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    family = family,
    controls_only = TRUE,
    v_folds = v_folds,
    method = Q_method
  )

  # re-define nodes to include RCT only data
  S <- data[, S_node]
  W <- data[S == 1, W_node]
  W <- cbind(W, procova_score)
  A <- data[S == 1, A_node]
  Y <- data[S == 1, Y_node]
  delta <- as.numeric(!is.na(Y))
  n <- nrow(W)

  # estimate ATE using TMLE
  Q <- learn_Q(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    method = Q_method,
    v_folds = 5,
    family = family,
    theta_bounds = c(-Inf, Inf)
  )

  # targeting Q
  Q_star <- tmle(
    Y = Y, A = A, W = W, g1W = rep(g_rct, sum(S == 1)),
    Q = as.matrix(data.frame(Q$A1, Q$A0)),
    Delta = delta,
    g.Delta.SL.library = c("SL.glm"),
    family = family
  )
  psi_est <- Q_star$estimates$ATE$psi
  psi_eic <- Q_star$estimates$IC$IC.ATE
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE) / n)
  psi_ci_lower <- psi_est - 1.96 * psi_se
  psi_ci_upper <- psi_est + 1.96 * psi_se

  return(list(
    est = psi_est,
    lower = psi_ci_lower,
    upper = psi_ci_upper
  ))
}
