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
  # Initialize components
  W_comp <- NULL  # W-component of the EIC
  Pi_comp <- NULL # Pi-component of the EIC
  beta_comp <- NULL  # beta-component of the EIC

  # Replace missing Y's with 0
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  # Define a function to compute the SVD-based pseudoinverse
  svd_pseudoinverse <- function(mat, tol = 1e-3) {
    svd_res <- svd(mat)
    # Invert only singular values above the tolerance; set the rest to zero
    D_inv <- ifelse(svd_res$d > tol, 1 / svd_res$d, 0)
    # Reconstruct the pseudoinverse
    return(svd_res$v %*% diag(D_inv) %*% t(svd_res$u))
  }

  if (controls_only) {
    W_comp <- (1 - Pi$A0) * tau$A0 - psi_pound_est
    Pi_comp <- -1 / (1 - g) * tau$A0 * (S - Pi$pred)

    # Compute the information matrix
    IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
    # Compute its pseudoinverse via SVD
    IM_inv <- svd_pseudoinverse(IM, tol = 1e-3)

    # Multiply the pseudoinverse with the column means of the auxiliary basis (weighted by (1-Pi$A0))
    IM_A0 <- IM_inv %*% colMeans(tau$x_basis_A0 * (1 - Pi$A0))
    beta_comp <- as.numeric(tau$x_basis %*% IM_A0) *
      (S - Pi$pred) *
      (Y_tmp - theta - (S - Pi$pred) * tau$A0) *
      weights
  } else {
    W_comp <- (1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1 - psi_pound_est
    Pi_comp <- (A / g * tau$A1 - (1 - A) / (1 - g) * tau$A0) * (S - Pi$pred)

    # Compute the information matrix
    IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
    # Use SVD to compute the pseudoinverse
    IM_inv <- svd_pseudoinverse(IM, tol = 1e-3)

    # Compute the D matrix used in the beta component
    D_mat <- tau$x_basis %*% IM_inv *
      (S - Pi$pred) *
      (Y_tmp - theta - (S - Pi$pred) * tau$pred) *
      weights

    # Compute beta component using the auxiliary bases
    if (ncol(D_mat) > 1) {
      beta_comp <- (rowSums(D_mat %*% diag(colMeans((1 - Pi$A0) * tau$x_basis_A0))) -
                      rowSums(D_mat %*% diag(colMeans((1 - Pi$A1) * tau$x_basis_A1))))
    } else {
      beta_comp <- (rowSums(D_mat * colMeans((1 - Pi$A0) * tau$x_basis_A0)) -
                      rowSums(D_mat * colMeans((1 - Pi$A1) * tau$x_basis_A1)))
    }
  }

  # Return the sum of all components as the efficient influence curve
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
# IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
# D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)

get_eic_psi_tilde <- function(psi_tilde, g, theta, Y, A, n, weights) {

  # Replace NA's in Y with 0
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  # Compute the information matrix using the original basis
  IM <- t(psi_tilde$x_basis) %*% diag(g * (1 - g)) %*% psi_tilde$x_basis / n

  # Use SVD to compute the pseudoinverse of IM
  svd_IM <- svd(IM)
  U <- svd_IM$u
  D <- svd_IM$d
  V <- svd_IM$v

  # Set a tolerance for small singular values
  tol <- 1e-3
  # Invert singular values above the tolerance; for values below, set the inverse to zero.
  D_inv <- ifelse(D > tol, 1 / D, 0)
  # Construct the pseudoinverse of IM using the SVD components
  IM_inv <- V %*% diag(D_inv) %*% t(U)

  # Compute the adjustment (beta) component of the efficient influence curve
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
