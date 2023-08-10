max_iter=10
eps=1e-5
tmle=FALSE
psi_tilde_adaptive=TRUE

atmle <- function(max_iter=10,
                  eps=1e-5,
                  external_trx=TRUE,
                  tmle=FALSE,
                  psi_tilde_adaptive=TRUE) {
  n <- length(A)

  # estimate bias psi_pound ----------------------------------------------------
  # learn relevant parts
  theta_pred <- learn_theta(W, A, Y)
  Pi_pred <- learn_Pi(S, W, A)
  tau_pred <- learn_tau(S, W, A, Y, Pi_pred$pred, theta_pred)
  g_pred <- learn_g(W, A)

  psi_pound_est <- NULL
  psi_pound_eic <- NULL

  if (tmle) {
    # perform TMLE update of Pi
    cur_iter <- 0
    tol <- Inf

    while (cur_iter <= max_iter && tol > eps) {
      print(cur_iter)
      # update
      Pi_star <- Pi_tmle(S, W, A, tau_pred, Pi_pred)

      # re-learn tau using Pi_star
      tau_star <- learn_tau(S, W, A, Y, Pi_star$pred, theta_pred)

      # current EIC
      cur_eic <- get_eic_Pi(S, A, g_pred, tau_star, Pi_star)
      tol <- abs(mean(cur_eic) - 0)

      cur_iter <- cur_iter + 1
      print(tol)
    }

    psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0)-mean((1-Pi_star$A1)*tau_star$A1)
    psi_pound_eic <- get_eic_psi_pound()

  } else {
    # no targeting, just plug-in
    psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0)-mean((1-Pi_pred$A1)*tau_pred$A1)
    psi_pound_eic <- get_eic_psi_pound()
  }

  # estimate pooled ATE psi_tilde ----------------------------------------------
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL

  if (psi_tilde_adaptive) {
    # use A-TMLE to estimate psi_tilde
    g_tilde <- learn_g(W, A)
    theta_tilde <- learn_theta(W, A, Y)
    psi_tilde <- learn_psi_tilde(W, A, Y, g_tilde, theta_tilde)

    psi_tilde_est <- mean(psi_tilde)
    psi_tilde_eic <- get_eic_psi_tilde()

  } else {
    # TODO: use TMLE to estimate psi_tilde
  }

  # estimate psi ---------------------------------------------------------------
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_pound_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci <- c(psi_est-1.96*psi_pound_se, psi_est+1.96*psi_pound_se)
}
