get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred)-mean(Pi$pred))
}

get_eic_psi_pound <- function(Pi, tau, g, theta, psi_pound_est, S, A, Y, n, controls_only) {
  if (controls_only) {
    W_comp <- (1-Pi$A0)*tau$A0-psi_pound_est
    Pi_comp <- -1/(1-g)*tau$A0*(S-Pi$pred)
    IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
    IM_A0 <- IM%*%colMeans(tau$x_basis_A0*(1-Pi$A0))
    tmp <- tau$A0
    beta_comp <- as.numeric(tau$x_basis%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)

    return(W_comp+Pi_comp+beta_comp)

  } else {
    W_comp <- (1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1-psi_pound_est
    Pi_comp <- ((A/g*tau$A1-(1-A)/(1-g)*tau$A0))*(S-Pi$pred)
    IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
    IM_A0 <- IM%*%colMeans(tau$x_basis_A0*(1-Pi$A0))
    IM_A1 <- IM%*%colMeans(tau$x_basis_A1*(1-Pi$A1))
    tmp <- vector(length = n)
    tmp[A == 1] <- tau$A1[A == 1]
    tmp[A == 0] <- tau$A0[A == 0]
    D_beta_A0 <- as.numeric(tau$x_basis%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
    D_beta_A1 <- as.numeric(tau$x_basis%*%IM_A1)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tmp)
    beta_comp <- D_beta_A0-D_beta_A1

    return(W_comp+Pi_comp+beta_comp)
  }
}

get_eic_psi_tilde <- function(psi_tilde, g_pred, theta, Y, A, n) {
  IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
  D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)
  return(as.vector(psi_tilde$pred-mean(psi_tilde$pred)+D_beta))
}

get_eic_psi_tilde_2 <- function(g, Q_A1, Q_A0, A, Y) {
  D_A1 <- A/g*(Y-Q_A1)
  D_A0 <- (1-A)/(1-g)*(Y-Q_A0)

  return(D_A1-D_A0+Q_A1-Q_A0-mean(Q_A1-Q_A0))
}

get_eic_psi_nonparametric <- function(Q, Pi, g, S, A, Y, psi_est) {
  W_comp <- Q$S1A1-Q$S1A0-psi_est
  Q_comp <- (S/Pi$pred)*(A/g-(1-A)/(1-g))*(Y-Q$pred)
  return(W_comp+Q_comp)
}
