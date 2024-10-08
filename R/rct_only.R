#' @title RCT only TMLE
#'
#' @export
rct_only <- function(data,
                     S,
                     W,
                     A,
                     Y,
                     g_rct,
                     family,
                     nuisance_method = "glm",
                     verbose = TRUE) {

  # define nodes
  S <- data[, S]
  W <- data[S == 1, W]
  A <- data[S == 1, A]
  Y <- data[S == 1, Y]
  delta <- as.numeric(!is.na(Y)) # missingness indicator
  n <- nrow(W)

  # estimate ATE using TMLE
  Q <- learn_Q(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    method = nuisance_method,
    v_folds = 5,
    family = family,
    theta_bounds = c(-Inf, Inf)
  )
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
