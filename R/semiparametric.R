library(hal9001)
library(data.table)

# SEMIPARAMETRIC ---------------------------------------------------------------
# LOSS FUNCTION DIRECTLY LEARN TAU ---------------------------------------------
# RELAXED HAL ------------------------------------------------------------------

# Psi_pound --------------------------------------------------------------------
# learn theta(W,A)=E(Y|W,A)
fit_theta <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), Y, "gaussian")
theta_pred <- fit_theta$pred

# learn Pi(1|W,A)=P(S=1|W,A)
fit_Pi <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), S, "binomial")
Pi_pred <- fit_Pi$pred

weights <- 1

# make pseudo outcome and weights
pseudo_outcome <- ifelse(abs(S - Pi_pred) < 1e-10, 0, (Y - theta_pred) / (S - Pi_pred))
pseudo_weights <- (S - Pi_pred)^2 * weights

# learn tau
keep <- which(abs(S - Pi_pred) > 1e-10)
fit_tau <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), pseudo_outcome[keep], "gaussian", weights = pseudo_weights[keep])
tau_pred <- fit_tau$pred
x_basis_A1 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 1)))
x_basis_A0 <- make_counter_design_matrix(fit_tau$basis_list, as.matrix(data.table(W, A = 0)))
tau_A1_pred <- x_basis_A1 %*% fit_tau$beta
tau_A0_pred <- x_basis_A0 %*% fit_tau$beta

# Psi_tilde --------------------------------------------------------------------

# learn theta_tilde(W)=E(Y|W)
fit_theta_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), Y, "gaussian")
theta_tilde_pred <- fit_theta_tilde$pred

# learn g(1|W)=P(A=1|W)
fit_g <- fit_relaxed_hal(as.matrix(data.table(W)), A, "binomial")
g_pred <- fit_g$pred

weights <- 1

# make pseudo outcome and weights
pseudo_outcome_tilde <- ifelse(abs(A - g_pred) < 1e-10, 0, (Y - theta_tilde_pred) / (A - g_pred))
pseudo_weights_tilde <- (A - g_pred)^2 * weights

# learn Psi_tilde
keep <- which(abs(A - g_pred) > 1e-10)
fit_psi_tilde <- fit_relaxed_hal(as.matrix(data.table(W)), pseudo_outcome_tilde[keep], "gaussian", weights = pseudo_weights_tilde[keep])
psi_tilde_pred <- mean(fit_psi_tilde$pred)

# target parameter point estimate
x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
Pi_A1_pred <- x_basis_A1 %*% fit_Pi$beta
Pi_A0_pred <- x_basis_A0 %*% fit_Pi$beta
psi_pound_pred <- mean((1 - Pi_A0_pred) * tau_A0_pred - (1 - Pi_A1_pred) * tau_A1_pred)
psi_pred <- psi_tilde_pred - psi_pound_pred
print(psi_pred)

