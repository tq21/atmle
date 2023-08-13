library(hal9001)
library(data.table)

# function to learn theta(W,A)=E(Y|W,A), relaxed HAL
learn_theta <- function(W, A, Y) {
  fit_theta <- fit_relaxed_hal(as.matrix(data.frame(W, A = A)), Y, "gaussian")
  pred <- as.vector(fit_theta$pred)

  return(pred)
}

# function to learn theta_tilde(W)=E(Y|W), relaxed HAL
learn_theta_tilde <- function(W, Y) {
  fit_theta <- fit_relaxed_hal(as.matrix(data.frame(W)), Y, "gaussian")
  pred <- as.vector(fit_theta$pred)

  return(pred)
}

# function to learn Pi(1|W,A)=P(S=1|W,A), relaxed HAL
learn_Pi <- function(S, W, A) {
  fit_Pi <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), S, "binomial")
  pred <- as.vector(fit_Pi$pred)
  x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- plogis(as.vector(x_basis_A1 %*% fit_Pi$beta))
  A0 <- plogis(as.vector(x_basis_A0 %*% fit_Pi$beta))

  return(list(pred = as.vector(pred),
              A1 = A1,
              A0 = A0))
}

# function to learn tau, R-learner, relaxed HAL
learn_tau <- function(S, W, A, Y, Pi, theta) {
  weights <- 1
  pseudo_outcome <- ifelse(abs(S-Pi) < 1e-3, 0, (Y-theta)/(S-Pi))
  pseudo_weights <- (S-Pi)^2*weights
  keep <- which(abs(S-Pi) > 1e-3)

  fit_tau <- fit_relaxed_hal(as.matrix(data.table(W, A = A)[keep, , drop = F]),
                             pseudo_outcome[keep], "gaussian", weights = pseudo_weights[keep])
  pred <- as.vector(fit_tau$pred)

  x_basis <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A)))
  x_basis_A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- as.vector(x_basis_A1 %*% fit_tau$beta)
  A0 <- as.vector(x_basis_A0 %*% fit_tau$beta)

  return(list(pred = pred,
              A1 = A1,
              A0 = A0,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0))
}

# function to learn psi_tilde, R-learner, relaxed HAL
learn_psi_tilde <- function(W, A, Y, g, theta) {
  weights <- 1
  pseudo_outcome <- (Y-theta)/(A-g)
  pseudo_weights <- (A-g)^2*weights

  fit_psi_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome, "gaussian", weights = pseudo_weights)
  pred <- as.vector(fit_psi_tilde$pred)

  x_basis <- make_counter_design_matrix(fit_psi_tilde$basis_list, as.matrix(data.table(W)))

  return(list(pred = pred,
              x_basis = x_basis))
}

# TODO: function to learn psi_tilde using TMLE
learn_psi_tilde_tmle <- function(g) {
  # clever covariates
  H1_n <- A/g*tau_pred$A1
  H0_n <- (1-A)/(1-g_pred)*tau_pred$A0

  # logistic submodel
  submodel <- glm(S ~ -1 + H1_n + H0_n + offset(Pi_pred$A1), family = "binomial")
  epsilon <- coef(submodel)

  # TMLE updates
  Pi_A1_star <- plogis(Pi_pred$A1 + epsilon["H1_n"] * H1_n)
  Pi_A0_star <- plogis(Pi_pred$A0 + epsilon["H0_n"] * H0_n)
  pred <- A * Pi_A1_star + (1 - A) * Pi_A0_star

  return(list(pred = as.vector(pred),
              A1 = as.vector(Pi_A1_star),
              A0 = as.vector(Pi_A0_star)))
}

# function to learn g(1|W)=P(1|W)
learn_g <- function(S, W, A, g_rct) {
  fit_g <- fit_relaxed_hal(as.matrix(data.table(W[S == 0,])), A[S == 0], "binomial",
                           smoothness_orders = 0, num_knots = 50)
  pred <- vector(length = length(S))
  pred[S == 1] <- g_rct
  pred[S == 0] <- as.vector(fit_g$pred)

  return(pred)
}

# function to perform TMLE update of Pi
Pi_tmle <- function(S, W, A, g, tau, Pi, target_gwt=FALSE) {

  # whether to target using weight
  if(target_gwt){
    wt <- A/g*tau$A1+(1-A)/(1-g)*tau$A0
    H1_n <- A
    H0_n <- 1-A
  } else{
    wt <- rep(1, length(A))
    H1_n <- A/g*tau$A1
    H0_n <- (1-A)/(1-g)*tau$A0
  }

  # logistic submodel
  suppressWarnings(
    epsilon <- coef(glm(S ~ -1 + offset(Pi_pred$A1) + H1_n + H0_n,
                        family = "binomial", weights = wt))
  )
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Pi_star <- NULL
  if (target_gwt) {
    Pi_star$pred <- plogis(Pi$pred+epsilon[1]*H1_n+epsilon[2]*H0_n)
    Pi_star$A1 <- plogis(Pi$pred+epsilon[1])
    Pi_star$A0 <- plogis(Pi$pred+epsilon[2])
  } else {
    Pi_star$pred <- plogis(Pi$pred+epsilon[1]*H1_n+epsilon[2]*H0_n)
    Pi_star$A1 <- plogis(Pi$pred+epsilon[1]/g)
    Pi_star$A0 <- plogis(Pi$pred+epsilon[2]/(1-g))
  }

  return(Pi_star)
}

get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred)-mean(Pi$pred))
}

get_eic_psi_pound_parametric <- function(Pi, tau, g, psi_pound, S, A, Y, n) {
  W_comp <- as.vector((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
  Pi_comp <- as.vector(A/g*tau$A1*(S-Pi$A1)-(1-A)/(1-g)*tau$A0*(S-Pi$A0))
  D_beta <- solve(t(tau$x_basis)%*%tau$x_basis/n)%*%t(tau$x_basis)%*%diag(Y-tau$pred)
  beta_comp <- as.vector((1-Pi$A0)%*%((tau$x_basis_S1A0-tau$x_basis_S0A0)/n)%*%D_beta-
    (1-Pi$A1)%*%((tau$x_basis_S1A1-tau$x_basis_S0A1)/n)%*%D_beta)
  return(W_comp+Pi_comp+beta_comp-psi_pound)
}


get_eic_psi_pound <- function(Pi, tau, g, theta, psi_pound, S, A, Y, n) {
  # TODO: get the EIC of psi_pound
  W_comp <- Pi$A0*tau$A0-Pi$A1*tau$A1
  Pi_comp <- (A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred)
  IM <- solve(t(tau$x_basis)%*%diag((Pi$pred*(1-Pi$pred)))%*%tau$x_basis/n)
  D_beta_A0 <- tau$x_basis%*%(IM%*%colMeans(tau$x_basis_A0))*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$pred)
  D_beta_A1 <- tau$x_basis%*%(IM%*%colMeans(tau$x_basis_A1))*(S-Pi$pred)*(Y-theta-(S-Pi$pred)*tau$pred)
  beta_comp <- D_beta_A0*Pi$A0-D_beta_A1*Pi$A1

  return(as.vector(W_comp+Pi_comp+beta_comp-psi_pound))
}

get_eic_psi_tilde <- function(psi_tilde, g_pred, theta, Y, A) {
  IM <- solve(t(psi_tilde$x_basis)%*%diag((g_pred*(1-g_pred)))%*%psi_tilde$x_basis/n)%*%colMeans(psi_tilde$x_basis)
  D_beta <- psi_tilde$x_basis%*%IM*(A-g_pred)*(Y-theta-(A-g_pred)*psi_tilde$pred)
  return(as.vector(psi_tilde$pred-mean(psi_tilde$pred)+D_beta))
}

learn_Q <- function(W, A, Y) {
  fit_Q <- fit_relaxed_hal(as.matrix(data.table(W, A)), Y, "gaussian")
  pred <- as.vector(fit_Q$pred)

  x_basis_A1 <- make_counter_design_matrix(fit_Q$basis_list, as.matrix(data.table(W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_Q$basis_list, as.matrix(data.table(W, A = 0)))
  A1 <- as.vector(x_basis_A1 %*% fit_Q$beta)
  A0 <- as.vector(x_basis_A0 %*% fit_Q$beta)

  return(list(pred = pred,
              A1 = A1,
              A0 = A0))
}

Q_tmle <- function(g, Q, A, Y_bound) {

  wt <- A/g+(1-A)/(1-g)
  H1W <- A
  H0W <- 1-A

  submodel <- glm(Y_bound ~ -1 + offset(Q$pred) + H0W + H1W, family = "quasibinomial", weights = wt)
  epsilon <- coef(submodel)

  Q_star <- Q$pred+epsilon[1]*H0W+epsilon[2]*H1W
  Q_A1_star <- Q$A1+rep(epsilon[1], length(Y_bound))
  Q_A0_star <- Q$A0+rep(epsilon[2], length(Y_bound))

  return(list(Q_star = Q_star,
              A1 = Q_A1_star,
              A0 = Q_A0_star))
}

get_eic_psi_tilde_2 <- function(g, Q_A1, Q_A0, A, Y) {
  D_A1 <- A/g*(Y-Q_A1)
  D_A0 <- (1-A)/(1-g)*(Y-Q_A0)

  return(D_A1-D_A0+Q_A1-Q_A0-mean(Q_A1-Q_A0))
}

bound <- function(X) {
  X_max <- max(X, na.rm = TRUE)
  X_min <- min(X, na.rm = TRUE)

  return((X-X_min)/(X_max-X_min))
}

learn_tau_parametric <- function(S, W, A, Y) {
  fit_tau <- fit_relaxed_hal(X = as.matrix(data.table(S, W, A)), Y = Y, family = "gaussian",
                             smoothness_orders = 1, num_knots = c(50,50,20))
  pred <- as.vector(fit_tau$pred)

  x_basis <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A)))
  x_basis_A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A = 1)))
  x_basis_A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S, W, A = 0)))
  x_basis_S1A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 1, W, A = 1)))
  x_basis_S0A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 0, W, A = 1)))
  x_basis_S1A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 1, W, A = 0)))
  x_basis_S0A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(S = 0, W, A = 0)))

  A1 <- as.vector(x_basis_S1A1 %*% fit_tau$beta) - as.vector(x_basis_S0A1 %*% fit_tau$beta)
  A0 <- as.vector(x_basis_S1A0 %*% fit_tau$beta) - as.vector(x_basis_S0A0 %*% fit_tau$beta)

  return(list(pred = pred,
              x_basis = x_basis,
              x_basis_A1 = x_basis_A1,
              x_basis_A0 = x_basis_A0,
              x_basis_S1A1 = x_basis_S1A1,
              x_basis_S0A1 = x_basis_S0A1,
              x_basis_S1A0 = x_basis_S1A0,
              x_basis_S0A0 = x_basis_S0A0,
              A1 = A1,
              A0 = A0))
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
