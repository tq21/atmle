max_iter=10
eps=1e-5

atmle <- function(max_iter=10,
                  eps=1e-5) {
  # estimate bias psi_pound ----------------------------------------------------
  # learn relevant parts
  theta_pred <- learn_theta(W, A, Y)
  Pi_pred <- learn_Pi(S, W, A)
  tau_pred <- learn_tau(S, W, A, Y, Pi_pred$pred, theta_pred)
  g_pred <- learn_g(W, A)

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
    cur_EIC_A1 <- get_EIC_Pi(S, A, g_pred, tau_star, Pi_star, A1 = TRUE)
    cur_EIC_A0 <- get_EIC_Pi(S, A, g_pred, tau_star, Pi_star, A1 = FALSE)
    tol <- abs(mean(c(cur_EIC_A1, cur_EIC_A0)) - 0)

    cur_iter <- cur_iter + 1
    print(tol)
  }

  EIC <- cur_EIC_A1 - cur_EIC_A0

  psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0)-mean((1-Pi_star$A1)*tau_star$A1)
  psi_pound_var <- var(EIC, na.rm = TRUE)/n

  # estimate pooled ATE psi_tilde
}







