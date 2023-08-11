max_iter=100
eps=1e-5
tmle=TRUE
psi_tilde_adaptive=TRUE

`%+%` <- function(a, b) paste0(a, b)

atmle <- function(max_iter=10,
                  eps=1e-5,
                  external_trx=TRUE,
                  tmle=FALSE,
                  psi_tilde_adaptive=TRUE) {
  n <- length(A)

  max_Y <- max(Y, na.rm = TRUE)
  min_Y <- min(Y, na.rm = TRUE)

  Y_bound <- bound(Y)

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
      cur_eic <- get_eic_Pi(g_pred, tau_star, Pi_star, S, A)
      tol <- abs(mean(cur_eic) - 0)

      cur_iter <- cur_iter + 1
      print(tol)
    }

    psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0-(1-Pi_star$A1)*tau_star$A1)
    psi_pound_eic <- get_eic_psi_pound(Pi_star, tau_star, g_pred, theta_pred, psi_pound_est, S, A, Y, n)

  } else {
    # no targeting, just plug-in
    psi_pound_est <- mean((1-Pi_pred$A0)*tau_pred$A0-(1-Pi_pred$A1)*tau_pred$A1)
    psi_pound_eic <- get_eic_psi_pound(Pi_pred, tau_pred, psi_pound_est, g_pred, theta_pred, S, A, Y, n)
  }

  # estimate pooled ATE psi_tilde ----------------------------------------------
  psi_tilde_est <- NULL
  psi_tilde_eic <- NULL

  if (psi_tilde_adaptive) {
    # use A-TMLE to estimate psi_tilde
    psi_tilde <- learn_psi_tilde(W, A, Y, g_pred, theta_pred)

    psi_tilde_est <- mean(psi_tilde$pred)
    psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g_pred, theta_pred, Y, A)

  } else {
    # use TMLE to estimate psi_tilde
    Q_pred <- learn_Q(W, A, Y)

    Q_star <- Q_tmle(g_pred, Q_pred, A, Y_bound)
    Q_A1 <- mean((max_Y-min_Y)*Q_star$A1+min_Y, na.rm = TRUE)
    Q_A0 <- mean((max_Y-min_Y)*Q_star$A0+min_Y, na.rm = TRUE)

    psi_tilde_est <- Q_A1 - Q_A0
    psi_tilde_eic <- get_eic_psi_tilde_2(g_pred, Q_A1, Q_A0, A, Y)
  }

  # estimate psi ---------------------------------------------------------------
  psi_est <- psi_tilde_est - psi_pound_est
  psi_eic <- psi_tilde_eic - psi_pound_eic
  psi_pound_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci <- as.character(psi_est-1.96*psi_pound_se) %+% ", " %+% as.character(psi_est+1.96*psi_pound_se)

  print("point estimate: " %+% psi_est %+% ", 95% CI: " %+% psi_ci)
}



atmle(max_iter=0,eps=1e-5,external_trx=TRUE,tmle=FALSE,psi_tilde_adaptive=FALSE)
atmle(max_iter=0,eps=1e-5,external_trx=TRUE,tmle=TRUE,psi_tilde_adaptive=FALSE)
atmle(max_iter=0,eps=1e-5,external_trx=TRUE,tmle=TRUE,psi_tilde_adaptive=TRUE)


