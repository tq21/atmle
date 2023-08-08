library(hal9001)
library(data.table)

# function to learn theta(W,A)=E(Y|W,A), relaxed HAL
learn_theta <- function(W, A, Y) {
  fit_theta <- fit_relaxed_hal(as.matrix(data.frame(W, A = A)), Y, "gaussian")
  pred <- as.vector(fit_theta$pred)

  return(pred)
}

# function to learn Pi(1|W,A)=P(S=1|W,A), relaxed HAL
learn_Pi <- function(S, W, A) {
  fit_Pi <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), S, "binomial")
  pred <- as.vector(fit_Pi$pred)
  x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- as.vector(x_basis_A1 %*% fit_Pi$beta)
  A0 <- as.vector(x_basis_A0 %*% fit_Pi$beta)

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}

# function to learn tau, Lars' loss function, relaxed HAL
learn_tau <- function(S, W, A, Y, Pi_pred, theta_pred) {
  weights <- 1
  pseudo_outcome <- ifelse(abs(S - Pi_pred) < 1e-10, 0, (Y - theta_pred) / (S - Pi_pred))
  pseudo_weights <- (S - Pi_pred)^2 * weights
  keep <- which(abs(S - Pi_pred) > 1e-10)
  fit_tau <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), pseudo_outcome[keep], "gaussian", weights = pseudo_weights[keep])
  pred <- as.vector(fit_tau$pred)
  x_basis_A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- as.vector(x_basis_A1 %*% fit_tau$beta)
  A0 <- as.vector(x_basis_A0 %*% fit_tau$beta)

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}

# function to learn g(1|W)=P(1|W)
learn_g <- function(W, A) {
  fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
  pred <- as.vector(fit_g$pred)

  return(pred)
}

# function to perform TMLE update of Pi
Pi_tmle <- function(S, W, A, tau_pred, Pi_pred) {
  # learn g(1|W)=P(1|W)
  fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
  g_pred <- fit_g$pred

  # clever covariates
  H1_n <- as.numeric(A == 1)/g_pred*tau_pred$A1
  H0_n <- -as.numeric(A == 0)/(1-g_pred)*tau_pred$A0

  # logistic submodel
  submodel_A1 <- glm(S ~ -1 + offset(Pi_pred$A1) + H1_n, family = "binomial")
  submodel_A0 <- glm(S ~ -1 + offset(Pi_pred$A0) + H0_n, family = "binomial")

  # TMLE updates
  Pi_A1_star <- as.vector(predict(submodel_A1, type = "response"))
  Pi_A0_star <- as.vector(predict(submodel_A0, type = "response"))

  pred <- vector(length = length(Pi_A1_star))
  pred[A == 1] <- Pi_A1_star[A == 1]
  pred[A == 0] <- Pi_A0_star[A == 0]

  return(list(pred = pred,
              A1 = Pi_A1_star,
              A0 = Pi_A0_star))
}

get_EIC_Pi <- function(S, A, g_pred, tau_pred, Pi_pred, A1) {
  if (A1) {
    return(1/g_pred*tau_pred$A1*(S-Pi_pred$A1))
  } else {
    return(1/(1-g_pred)*tau_pred$A0*(S-Pi_pred$A0))
  }
}

# psi_pound_pred <- mean((1-Pi_star$Pi_A0_star)*tau_star$tau_A0_star)-mean((1-Pi_star$Pi_A1_star)*tau_star$tau_A1_star)
#

#
# # Psi_tilde --------------------------------------------------------------------
#
# # learn theta_tilde(W)=E(Y|W)
# fit_theta_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), Y, "gaussian")
# theta_tilde_pred <- fit_theta_tilde$pred
#
# # learn g(1|W)=P(A=1|W)
# fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
# g_pred <- fit_g$pred
#
# weights <- 1
#
# # make pseudo outcome and weights
# pseudo_outcome_tilde <- ifelse(abs(A - g_pred) < 1e-10, 0, (Y - theta_tilde_pred) / (A - g_pred))
# pseudo_weights_tilde <- (A - g_pred)^2 * weights
#
# # learn Psi_tilde
# keep <- which(abs(A - g_pred) > 1e-10)
# fit_psi_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome_tilde[keep], "gaussian", weights = pseudo_weights_tilde[keep])
# psi_tilde_pred <- mean(fit_psi_tilde$pred)
#
# # target parameter point estimate
# x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
# x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
# Pi_A1_pred <- x_basis_A1 %*% fit_Pi$beta
# Pi_A0_pred <- x_basis_A0 %*% fit_Pi$beta
# psi_pound_pred <- mean((1 - Pi_A0_pred) * tau_A0_pred - (1 - Pi_A1_pred) * tau_A1_pred)
# psi_pred <- psi_tilde_pred - psi_pound_pred
# print(psi_pred)
