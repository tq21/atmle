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
# get_eic_psi_pound <- function(Pi,
#                               tau,
#                               g,
#                               theta,
#                               psi_pound_est,
#                               S,
#                               A,
#                               Y,
#                               n,
#                               controls_only,
#                               weights,
#                               eic_method = "svd_pseudo_inv") {
#
#   Y_tmp <- Y
#   Y_tmp[is.na(Y)] <- 0
#
#   if (controls_only) {
#     W_comp <- (1 - Pi$A0) * tau$A0 - psi_pound_est
#     Pi_comp <- -1 / (1 - g) * tau$A0 * (S - Pi$pred)
#     IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
#     if (dim(tau$x_basis)[2] == 1) {
#       IM_inv <- solve(IM)
#     } else {
#       if (eic_method == "svd_pseudo_inv") {
#         # SVD-based pseudo-inverse
#         IM_inv <- svd_pseudo_inv(IM)
#       } else if (eic_method == "diag") {
#         IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
#       }
#     }
#     IM_A0 <- IM_inv %*% colMeans(tau$x_basis_A0 * (1 - Pi$A0))
#     beta_comp <- as.numeric(tau$x_basis %*% IM_A0) *
#       (S - Pi$pred) *
#       (Y_tmp - theta - (S - Pi$pred) * tau$A0) *
#       weights
#   } else {
#     W_comp <- (1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1 - psi_pound_est
#     Pi_comp <- (A / g * tau$A1 - (1 - A) / (1 - g) * tau$A0) * (S - Pi$pred)
#     IM <- t(tau$x_basis) %*% diag(Pi$pred * (1 - Pi$pred)) %*% tau$x_basis / n
#     if (dim(tau$x_basis)[2] == 1) {
#       IM_inv <- solve(IM)
#     } else {
#       if (eic_method == "svd_pseudo_inv") {
#         # SVD-based pseudo-inverse
#         IM_inv <- svd_pseudo_inv(IM)
#       } else if (eic_method == "diag") {
#         IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
#       }
#     }
#     D_mat <- tau$x_basis %*% IM_inv *
#       (S - Pi$pred) *
#       (Y_tmp - theta - (S - Pi$pred) * tau$pred) *
#       weights
#     if (ncol(D_mat) > 1) {
#       beta_comp <- (rowSums(D_mat %*% diag(colMeans((1 - Pi$A0) * tau$x_basis_A0))) -
#                       rowSums(D_mat %*% diag(colMeans((1 - Pi$A1) * tau$x_basis_A1))))
#     } else {
#       beta_comp <- (rowSums(D_mat * colMeans((1 - Pi$A0) * tau$x_basis_A0)) -
#                       rowSums(D_mat * colMeans((1 - Pi$A1) * tau$x_basis_A1)))
#     }
#   }
#
#   return(W_comp + Pi_comp + beta_comp)
# }

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
  W_comp <- NULL # W-component of the EIC
  Pi_comp <- NULL # Pi-component of the EIC
  beta_comp <- NULL # beta-component of the EIC

  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  if (controls_only) {
    W_comp <- (1 - Pi$A0) * tau$A0 - psi_pound_est
    Pi_comp <- -1 / (1 - g) * tau$A0 * (S - Pi$pred)
    IM <- solve(t(tau$x_basis) %*% diag((Pi$pred * (1 - Pi$pred))) %*% tau$x_basis / n)
    IM_A0 <- IM %*% colMeans(tau$x_basis_A0 * (1 - Pi$A0))
    beta_comp <- as.numeric(tau$x_basis %*% IM_A0) * (S - Pi$pred) * (Y_tmp - theta - (S - Pi$pred) * tau$A0) * weights
  } else {
    W_comp <- (1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1 - psi_pound_est
    Pi_comp <- (A / g * tau$A1 - (1 - A) / (1 - g) * tau$A0) * (S - Pi$pred)
    IM <- t(tau$x_basis) %*% diag((Pi$pred * (1 - Pi$pred))) %*% tau$x_basis / n
    D <- tau$x_basis %*% solve(IM) * (S - Pi$pred) * (Y_tmp - theta - (S - Pi$pred) * tau$pred) * weights
    beta_comp <- NULL
    if (ncol(D) > 1) {
      beta_comp <- (rowSums(D %*% diag(colMeans((1 - Pi$A0) * tau$x_basis_A0))) - rowSums(D %*% diag(colMeans((1 - Pi$A1) * tau$x_basis_A1))))
    } else {
      beta_comp <- (rowSums(D * colMeans((1 - Pi$A0) * tau$x_basis_A0)) - rowSums(D * colMeans((1 - Pi$A1) * tau$x_basis_A1)))
    }
  }
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
get_eic_psi_tilde <- function(psi_tilde, g, theta, Y, A, n, weights, eic_method = "svd_pseudo_inv") {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0
  IM <- t(psi_tilde$x_basis) %*% diag(g * (1 - g)) %*% psi_tilde$x_basis / n
  if (dim(psi_tilde$x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
    }
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

get_beta_h_T <- function(x_basis,
                         g1W,
                         eic_method = "svd_pseudo_inv") {
  n <- nrow(x_basis)
  IM <- t(x_basis) %*% diag(g1W*(1-g1W)) %*% x_basis / n
  if (dim(x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
    }
  }
  beta_h <- as.vector(colMeans(x_basis) %*% IM_inv)
  return(beta_h)
}

get_beta_h <- function(x_basis,
                       x_basis_A0,
                       x_basis_A1,
                       Pi,
                       eic_method = "svd_pseudo_inv") {

  n <- nrow(x_basis)
  IM <- t(x_basis) %*% diag((Pi$pred * (1 - Pi$pred))) %*% x_basis / n
  if (dim(x_basis)[2] == 1) {
    IM_inv <- solve(IM)
  } else {
    if (eic_method == "svd_pseudo_inv") {
      # SVD-based pseudo-inverse
      IM_inv <- svd_pseudo_inv(IM)
    } else if (eic_method == "diag") {
      IM_inv <- solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
    }
  }
  beta_h <- as.vector(colMeans(x_basis_A0*(1-Pi$A0)-x_basis_A1*(1-Pi$A1)) %*% IM_inv)
  return(beta_h)
}

eic_psi_pound_wm <- function(S,
                             Y,
                             A,
                             g1W,
                             theta_WA,
                             Pi,
                             tau_S,
                             weights,
                             controls_only,
                             std_wrt_rct_W,
                             IM_inv = NULL,
                             eic_method = "svd_pseudo_inv") {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  if (is.null(IM_inv)) {
    IM <- t(tau_S$phi_WA) %*% diag(Pi$pred*(1-Pi$pred)) %*% tau_S$phi_WA / length(Y)
    IM_inv <- tryCatch({
      solve(IM)
    }, error = function(e) {
      if (eic_method == "svd_pseudo_inv") {
        svd_pseudo_inv(IM)
      } else if (eic_method == "diag") {
        solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
      } else {
        stop("Unknown eic_method specified.")
      }
    })
  }

  if (controls_only) {
    if (std_wrt_rct_W) {
      psi_pound_est <- mean((S/mean(S))*(1-Pi$A0)*tau_S$cate_W0)
      W_comp <- (S/mean(S))*((1-Pi$A0)*tau_S$cate_W0)-psi_pound_est
      Pi_comp <- (-(1-A)/(1-g1W)*tau_S$cate_W0)*(S-Pi$pred)
      D <- tau_S$phi_WA %*% IM_inv*(S-Pi$pred)*(Y_tmp-theta_WA-(S-Pi$pred)*tau_S$cate_WA)*weights
      if (ncol(D) > 1) {
        beta_comp <- (rowSums(D %*% diag(colMeans((S/mean(S))*(1-Pi$A0)*tau_S$phi_W0))))
      } else {
        beta_comp <- (rowSums(D * colMeans((S/mean(S))*(1-Pi$A0)*tau_S$phi_W0)))
      }
    } else {
      psi_pound_est <- mean((1-Pi$A0)*tau_S$cate_W0)
      W_comp <- (1-Pi$A0)*tau_S$cate_W0-psi_pound_est
      Pi_comp <- -(1-A)/(1-g1W)*tau_S$cate_W0*(S-Pi$pred)
      D <- tau_S$phi_WA %*% IM_inv*(S-Pi$pred)*(Y_tmp-theta_WA-(S-Pi$pred)*tau_S$cate_WA)*weights
      if (ncol(D) > 1) {
        beta_comp <- (rowSums(D %*% diag(colMeans((1-Pi$A0)*tau_S$phi_W0))))
      } else {
        beta_comp <- (rowSums(D * colMeans((1-Pi$A0)*tau_S$phi_W0)))
      }
    }
  } else {
    if (std_wrt_rct_W) {
      psi_pound_est <- mean((S/mean(S))*((1-Pi$A0)*tau_S$cate_W0-(1-Pi$A1)*tau_S$cate_W1))
      W_comp <- (S/mean(S))*((1-Pi$A0)*tau_S$cate_W0-(1-Pi$A1)*tau_S$cate_W1-psi_pound_est)
      Pi_comp <- (A/g1W*tau_S$cate_W1-(1-A)/(1-g1W)*tau_S$cate_W0)*(S-Pi$pred)
      D <- tau_S$phi_WA %*% IM_inv*(S-Pi$pred)*(Y_tmp-theta_WA-(S-Pi$pred)*tau_S$cate_WA)*weights
      if (ncol(D) > 1) {
        beta_comp <- (rowSums(D %*% diag(colMeans((S/mean(S))*(1-Pi$A0)*tau_S$phi_W0)))-rowSums(D %*% diag(colMeans((S/mean(S))*(1-Pi$A1)*tau_S$phi_W1))))
      } else {
        beta_comp <- (rowSums(D * colMeans((S/mean(S))*(1-Pi$A0)*tau_S$phi_W0))-rowSums(D * colMeans((S/mean(S))*(1-Pi$A1)*tau_S$phi_W1)))
      }
    } else {
      psi_pound_est <- mean((1-Pi$A0)*tau_S$cate_W0-(1-Pi$A1)*tau_S$cate_W1)
      W_comp <- (1-Pi$A0)*tau_S$cate_W0-(1-Pi$A1)*tau_S$cate_W1-psi_pound_est
      Pi_comp <- (A/g1W*tau_S$cate_W1-(1-A)/(1-g1W)*tau_S$cate_W0)*(S-Pi$pred)
      D <- tau_S$phi_WA %*% IM_inv*(S-Pi$pred)*(Y_tmp-theta_WA-(S-Pi$pred)*tau_S$cate_WA)*weights
      if (ncol(D) > 1) {
        beta_comp <- (rowSums(D %*% diag(colMeans((1-Pi$A0)*tau_S$phi_W0)))-rowSums(D %*% diag(colMeans((1-Pi$A1)*tau_S$phi_W1))))
      } else {
        beta_comp <- (rowSums(D * colMeans((1-Pi$A0)*tau_S$phi_W0))-rowSums(D * colMeans((1-Pi$A1)*tau_S$phi_W1)))
      }
    }
  }

  #print(paste("W_comp", mean(W_comp)))
  #print(paste("Pi_comp", mean(Pi_comp)))
  #print(paste("beta_comp", mean(beta_comp)))

  return(W_comp+Pi_comp+beta_comp)
}

eic_psi_tilde_wm <- function(S,
                             Y,
                             A,
                             g1W,
                             theta_W,
                             tau_A,
                             weights,
                             eic_method,
                             std_wrt_rct_W,
                             IM_inv = NULL) {
  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  if (is.null(IM_inv)) {
    IM <- t(tau_A$phi_W) %*% diag(g1W*(1-g1W)) %*% tau_A$phi_W / length(Y)
    IM_inv <- tryCatch({
      solve(IM)
    }, error = function(e) {
      if (eic_method == "svd_pseudo_inv") {
        svd_pseudo_inv(IM)
      } else if (eic_method == "diag") {
        solve(IM + diag(1e-3, nrow(IM), ncol(IM)))
      } else {
        stop("Unknown eic_method specified.")
      }
    })
  }
  if (std_wrt_rct_W) {
    D_beta <- weights*as.vector(tau_A$phi_W %*% IM_inv %*% colMeans((S/mean(S))*tau_A$phi_W) *
                                  (A-g1W)*(Y_tmp-theta_W-(A-g1W)*tau_A$cate))
    W_comp <- S/mean(S)*((tau_A$cate)-mean((S/mean(S))*tau_A$cate))
  } else {
    D_beta <- weights*as.vector(tau_A$phi_W %*% IM_inv %*% colMeans(tau_A$phi_W) *
                                  (A-g1W)*(Y_tmp-theta_W-(A-g1W)*tau_A$cate))
    W_comp <- tau_A$cate-mean(tau_A$cate)
  }

  #print(paste("W_comp", mean(W_comp)))
  #print(paste("D_beta", mean(D_beta)))

  return(W_comp+D_beta)
}
