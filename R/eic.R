#' @title Efficient influence function for the bias projection parameter
#' @description Computes the efficient influence function for the bias
#' projection parameter.
#'
#' @keywords internal
#'
#' @param Pi Results from \code{\link{learn_Pi}}.
#' @param tau Results from \code{\link{learn_tau}}.
#' @param g Results from \code{\link{learn_g}}.
#' @param theta Results from \code{\link{learn_theta}}.
#' @param psi_pound_est Estimate for the bias estimand.
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param Y A vector of outcomes.
#' @param n The number of observations.
#' @param controls_only A logical indicating whether external data has only
#' control-arm observations.
#' @param weights A vector of (e.g. inverse-censoring) weights.
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
                              controls_only,
                              weights) {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  if (controls_only) {
    W_comp <- (1 - Pi$A0) * tau$A0 - psi_pound_est
    Pi_comp <- -1 / (1 - g) * tau$A0 * (S - Pi$pred)
    IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
    if (dim(tau$x_basis)[2] == 1) {
      IM_inv <- solve(IM)
    } else {
      IM_inv <- svd_pseudo_inv(IM)
    }
    IM_A0 <- IM_inv %*% colMeans(tau$x_basis_A0 * (1 - Pi$A0))
    beta_comp <- as.numeric(tau$x_basis %*% IM_A0) *
      (S - Pi$pred) *
      (Y_tmp - theta - (S - Pi$pred) * tau$A0) *
      weights
  } else {
    W_comp <- (1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1 - psi_pound_est
    Pi_comp <- (A / g * tau$A1 - (1 - A) / (1 - g) * tau$A0) * (S - Pi$pred)
    IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
    if (dim(tau$x_basis)[2] == 1) {
      IM_inv <- solve(IM)
    } else {
      IM_inv <- svd_pseudo_inv(IM)
    }
    D_mat <- tau$x_basis %*% IM_inv *
      (S - Pi$pred) *
      (Y_tmp - theta - (S - Pi$pred) * tau$pred) *
      weights
    if (ncol(D_mat) > 1) {
      beta_comp <- (rowSums(D_mat %*% diag(colMeans((1 - Pi$A0) * tau$x_basis_A0))) -
                      rowSums(D_mat %*% diag(colMeans((1 - Pi$A1) * tau$x_basis_A1))))
    } else {
      beta_comp <- (rowSums(D_mat * colMeans((1 - Pi$A0) * tau$x_basis_A0)) -
                      rowSums(D_mat * colMeans((1 - Pi$A1) * tau$x_basis_A1)))
    }
  }

  return(W_comp + Pi_comp + beta_comp)
}

#' @title Efficient influence function for the pooled-ATE projection parameter
#' @description Computes the efficient influence function for the pooled-ATE
#' projection parameter.
#'
#' @keywords internal
#'
#' @param psi_tilde Results from \code{\link{learn_psi_tilde}}.
#' @param g_pred Results from \code{\link{learn_g}}.
#' @param theta Results from \code{\link{learn_theta}}.
#' @param Y A vector of outcomes.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param n The number of observations.
#' @param weights A vector of (e.g. inverse-censoring) weights.
#'
#' @return A vector of efficient influence function values.
get_eic_psi_tilde <- function(psi_tilde, g, theta, Y, A, n, weights) {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0
  IM <- t(psi_tilde$x_basis) %*% diag(g * (1 - g)) %*% psi_tilde$x_basis / n
  if (dim(psi_tilde$x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    # SVD-based pseudo-inverse
    IM_inv <- svd_pseudo_inv(IM)
  }
  D_beta <- weights * as.vector(
    psi_tilde$x_basis %*% IM_inv %*% colMeans(psi_tilde$x_basis) *
      (A - g) * (Y_tmp - theta - (A - g) * psi_tilde$pred)
  )
  W_comp <- psi_tilde$pred - mean(psi_tilde$pred)

  return(W_comp + D_beta)
}

get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A / g * tau$A1 - (1 - A) / (1 - g) * tau$A0) * (S - Pi$pred) - mean(Pi$pred))
}

get_eic_psi_tilde_2 <- function(g, Q_A1, Q_A0, A, Y) {
  D_A1 <- A / g * (Y - Q_A1)
  D_A0 <- (1 - A) / (1 - g) * (Y - Q_A0)

  return(D_A1 - D_A0 + Q_A1 - Q_A0 - mean(Q_A1 - Q_A0))
}

get_eic_psi_nonparametric <- function(Q, Pi, g, S, A, Y, psi_est, weights) {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0
  W_comp <- Q$S1A1 - Q$S1A0 - psi_est
  Q_comp <- (S / Pi$pred) * (A / g - (1 - A) / (1 - g)) * weights * (Y_tmp - Q$pred)
  return(W_comp + Q_comp)
}
