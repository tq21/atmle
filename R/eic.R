#' @title Efficient influence function for the bias parameter
#' @description Computes the efficient influence function for the bias parameter
#' under the working model.
#'
#' @keywords internal
#'
#' @param Pi Returns from \code{\link{learn_Pi}}.
#' @param tau Returns from \code{\link{learn_tau}}.
#' @param g Returns from \code{\link{learn_g}}.
#' @param theta Returns from \code{\link{learn_theta}}.
#' @param psi_pound_est Estimate of the bias estimand.
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param n The number of observations.
#' @param controls_only A logical indicating whether external data has only
#' control-arm observations.
#'
#' @return A vector of efficient influence function values.
get_eic_psi_pound <- function(Pi,
                              tau,
                              g,
                              theta,
                              psi_pound_est,
                              S,
                              A,
                              Y,
                              n,
                              controls_only) {

  W_comp <- NULL # W-component of the EIC
  Pi_comp <- NULL # Pi-component of the EIC
  beta_comp <- NULL # beta-component of the EIC

  if (controls_only) {
    W_comp <- (1-Pi$A0)*tau$A0-psi_pound_est
    Pi_comp <- -1/(1-g)*tau$A0*(S-Pi$pred)
    IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
    IM_A0 <- IM%*%colMeans(tau$x_basis_A0*(1-Pi$A0))
    beta_comp <- as.numeric(tau$x_basis%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$A0)

  } else {
    # W_comp <- (1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1-psi_pound_est
    # Pi_comp <- ((A/g*tau$A1-(1-A)/(1-g)*tau$A0))*(S-Pi$pred)
    # IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
    # IM_A0 <- IM%*%colMeans(tau$x_basis_A0*(1-Pi$A0))
    # IM_A1 <- IM%*%colMeans(tau$x_basis_A1*(1-Pi$A1))
    # D_beta_A0 <- as.numeric(tau$x_basis%*%IM_A0)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$pred)
    # D_beta_A1 <- as.numeric(tau$x_basis%*%IM_A1)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$pred)
    # beta_comp <- D_beta_A0-D_beta_A1

    #browser()

    W_comp <- (1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1-psi_pound_est
    Pi_comp <- ((A/g*tau$A1-(1-A)/(1-g)*tau$A0))*(S-Pi$pred)
    IM <- t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n
    D <- tau$x_basis%*%solve(IM)*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$pred)
    beta_comp <- rowSums(D*Pi$A0*tau$x_basis_A0)-rowSums(D*Pi$A1*tau$x_basis_A1)
  }

  # if (verbose) {
  #   print("EIC W component: " %+% round(mean(W_comp), 5))
  #   print("EIC Pi component: " %+% round(mean(Pi_comp), 5))
  #   print("EIC beta component: " %+% round(mean(beta_comp), 5))
  # }

  return(W_comp+Pi_comp+beta_comp)
}

#' @title Efficient influence function for the pooled-ATE parameter under the
#' working model.
#' @description Computes the efficient influence function for the pooled-ATE
#' parameter under the working model.
#'
#' @keywords internal
#'
#' @param psi_tilde Returns from \code{\link{learn_psi_tilde}}.
#' @param g_pred Returns from \code{\link{learn_g}}.
#' @param theta Returns from \code{\link{learn_theta}}.
#' @param Y A vector of outcomes.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param n The number of observations.
#'
#' @return A vector of efficient influence function values.
get_eic_psi_tilde <- function(psi_tilde,
                              g_pred,
                              theta,
                              Y,
                              A,
                              n) {
  # IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
  # D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)
  #browser()

  IM <- t(psi_tilde$x_basis)%*%diag(g_pred*(1-g_pred))%*%psi_tilde$x_basis/n
  D_beta <- psi_tilde$x_basis%*%solve(IM)*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)
  return(psi_tilde$pred-mean(psi_tilde$pred)+rowSums(D_beta*psi_tilde$x_basis))

  #return(as.vector(psi_tilde$pred-mean(psi_tilde$pred)+D_beta))
}

get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred)-mean(Pi$pred))
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
