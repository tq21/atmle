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

  # if (verbose) {
  # print("EIC W component: " %+% round(mean(W_comp), 5))
  # print("EIC Pi component: " %+% round(mean(Pi_comp), 5))
  # print("EIC beta component: " %+% round(mean(beta_comp), 5))
  # }

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
get_eic_psi_tilde <- function(psi_tilde,
                              g,
                              theta,
                              Y,
                              A,
                              n,
                              weights) {
  # IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
  # D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)

  Y_tmp <- Y
  Y_tmp[is.na(Y)] <- 0

  IM <- t(psi_tilde$x_basis) %*% diag(g * (1 - g)) %*% psi_tilde$x_basis / n
  D_beta <- weights * as.vector(psi_tilde$x_basis %*% solve(IM) %*% colMeans(psi_tilde$x_basis) * (A - g) * (Y_tmp - theta - (A - g) * psi_tilde$pred))
  W_comp <- psi_tilde$pred - mean(psi_tilde$pred)

  # print("Pooled EIC W component: " %+% round(mean(W_comp), 5))
  # print("Pooled EIC beta component: " %+% round(mean(D_beta), 5))

  return(W_comp + D_beta)

  # return(as.vector(psi_tilde$pred-mean(psi_tilde$pred)+D_beta))
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

get_eic_psi_tilde_surv <- function(data,
                                   data_long,
                                   g,
                                   Q_bar_r,
                                   stablize_weight,
                                   cate_surv,
                                   unique_t,
                                   Y,
                                   n) {

  # TODO: write this up

  tmp_1 <- stablize_weight*cate_surv$x_basis
  D_tilde_W_1 <- tmp_1*(Q_bar_r-cate_surv$pred)
  D_tilde_W_2 <- -1/n*tmp_1*Q_bar_r
  D_tilde_W <- D_tilde_W_1+D_tilde_W_2

  data_long[, D_tilde_Q := (clever_cov_A1-clever_cov_A0)*(dN_t-lambda)]
  D_tilde_Q <- map(unique_t, function(cur_t) {
    return(data_long[t == cur_t]$D_tilde_Q)
  })
  D_tilde_Q <- Reduce(`+`, D_tilde_Q)

  IM <- t(cate_surv$x_basis) %*% diag(g*(1-g)) %*% cate_surv$x_basis / n
  D_star_beta <- D_tilde_W %*% solve(IM) %*% colMeans(cate_surv$x_basis) + D_tilde_Q
  D_star <- cate_surv$pred - mean(cate_surv$pred) + as.numeric(D_star_beta)

  return(D_star)
}

get_eic_ipcw_r_learner <- function(data,
                                   Y,
                                   A,
                                   cate_surv,
                                   g,
                                   theta,
                                   n,
                                   weights) {

  IM <- t(cate_surv$x_basis)%*%diag(g*(1-g))%*%cate_surv$x_basis/n
  D_beta <- weights*as.vector(cate_surv$x_basis%*%solve(IM)%*%colMeans(cate_surv$x_basis)*(data[[A]]-g)*(data[[Y]]-theta-(data[[A]]-g)*cate_surv$pred))
  W_comp <- cate_surv$pred - mean(cate_surv$pred)

  return(W_comp + D_beta)
}

get_eic_surv_tmle <- function(data_long,
                              unique_t,
                              tmle) {

  data_long[, D_tilde_Q := (clever_cov_A1-clever_cov_A0)*(dN_t-lambda)]
  D_tilde_Q <- map(unique_t, function(cur_t) {
    return(data_long[t == cur_t]$D_tilde_Q)
  })
  D_tilde_Q <- Reduce(`+`, D_tilde_Q)
  D_star <- D_tilde_Q+data_long[t == t0]$surv_A1-data_long[t == t0]$surv_A0-tmle

  return(D_star)
}
