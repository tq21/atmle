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
  A1 <- plogis(as.vector(x_basis_A1 %*% fit_Pi$beta))
  A0 <- plogis(as.vector(x_basis_A0 %*% fit_Pi$beta))

  return(list(pred = as.vector(pred),
              A1 = A1,
              A0 = A0))
}

# function to learn tau, R-learner, relaxed HAL
learn_tau <- function(S, W, A, Y, Pi, theta) {
  weights <- 1
  pseudo_outcome <- (Y - theta) / (S - Pi)
  pseudo_weights <- (S - Pi)^2 * weights

  fit_tau <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), pseudo_outcome, "gaussian", weights = pseudo_weights)
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
  browser()

  weights <- 1
  pseudo_outcome <- ifelse(abs(A - g) < 1e-10, 0, (Y - theta) / (A - g))
  pseudo_weights <- (A - g)^2 * weights

  # fit hal
  hal_fit <- fit_hal(X = W, Y = pseudo_outcome, family = "gaussian", weights = pseudo_weights, smoothness_orders = 0)
  basis_list <- hal_fit$basis_list[hal_fit$coefs[-1] != 0]
  x_basis <- cbind(1, as.matrix(hal9001::make_design_matrix(W, basis_list)))

  # relaxed fit
  beta <- NULL
  pred <- NULL
  hal_relaxed_fit <- glm.fit(x = x_basis, y = Y, family = gaussian(), weights = weights, intercept = FALSE)
  beta <- coef(hal_relaxed_fit)
  beta[is.na(beta)] <- 0
  pred <- as.vector(x_basis %*% beta)


  coefs <- beta[beta != 0]
  tau <- x_basis %*% coef(glm.fit(x_basis, pseudo_outcome, weights = pseudo_weights))

  x_basis_proj <- model.matrix(~1, data = as.data.frame(W))
  coef_proj <-  coef(glm.fit(x_basis_proj, tau))
  tau_proj <- x_basis_proj %*% coef_proj

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
  fit_g <- fit_relaxed_hal(as.matrix(data.table(W[S == 0,])), A[S == 0], "binomial")
  pred <- c(rep(g_rct, sum(S == 1)), as.vector(fit_g$pred))

  return(pred)
}

# function to perform TMLE update of Pi
Pi_tmle <- function(S, W, A, tau_pred, Pi_pred) {
  # learn g(1|W)=P(1|W)
  fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
  g_pred <- fit_g$pred

  # clever covariates
  H1_n <- A/g_pred*tau_pred$A1
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

get_eic_Pi <- function(g, tau, Pi, S, A) {
  return((A/g*tau$A1-(1-A)/(1-g)*tau$A0)*(S-Pi$pred))
}

get_eic_psi_pound <- function(Pi, tau, g, theta, psi_pound, S, A, Y, n) {
  # TODO: get the EIC of psi_pound
  W_comp <- Pi$A0*tau$A0-Pi$A1*tau$A1
  Pi_comp <- get_eic_Pi(g, tau, Pi, S, A)
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

  browser()

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
