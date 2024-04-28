# A-TMLE for psi_tilde
psi_tilde_only_atmle <- function(data,
                                 S_node,
                                 W_node,
                                 A_node,
                                 Y_node,
                                 nuisance_method = "lasso",
                                 working_model = "lasso",
                                 p_rct = 0.5,
                                 verbose = TRUE,
                                 transform = TRUE) {

  # define nodes
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]
  n <- nrow(data)

  # estimate pooled ATE psi_tilde ----------------------------------------------
  # learn nuisance parts
  if (verbose) print("learning P(A|W)")
  g <- learn_g(S, W, A, p_rct, nuisance_method)

  if (verbose) print("learning E(Y|W)")
  theta_tilde <- learn_theta_tilde(W, Y, nuisance_method)

  # learn psi_tilde using R-loss
  if (verbose) print("learning psi_tilde")
  psi_tilde <- NULL
  if (transform) {
    psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, working_model)
  } else {
    psi_tilde <- learn_psi_tilde_test(W, A, Y, g, theta_tilde, working_model)
  }
  psi_tilde_est <- mean(psi_tilde$pred)
  psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g, theta_tilde, Y, A, n)

  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
  psi_tilde_ci_lower <- psi_tilde_est - 1.96 * psi_tilde_se
  psi_tilde_ci_upper <- psi_tilde_est + 1.96 * psi_tilde_se

  return(list(
    est = psi_tilde_est,
    lower = psi_tilde_ci_lower,
    upper = psi_tilde_ci_upper
  ))
}
