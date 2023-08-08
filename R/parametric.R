# library(hal9001)
# library(data.table)
#
# # PARAMETRIC -------------------------------------------------------------------
# # WORKING MODEL FOR TAU, SAME WORKING MODEL FOR PSI_TILDE ----------------------
# # RELAXED HAL ------------------------------------------------------------------
#
# # learn Pi(1|W,A)=P(S=1|W,A)
# fit_Pi <- fit_relaxed_hal(as.matrix(data.table(W, A = A)), S, "binomial")
# x_basis_A1 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 1)))
# x_basis_A0 <- make_counter_design_matrix(fit_Pi$basis_list, as.matrix(data.table(W, A = 0)))
# Pi_A1_pred <- x_basis_A1 %*% fit_Pi$beta
# Pi_A0_pred <- x_basis_A0 %*% fit_Pi$beta
#
# # learn Q_bar(S,W,A)=E(Y|S,W,A)
# fit_Q_bar <- fit_relaxed_hal(as.matrix(data.table(S = S, W, A = A)), Y, "gaussian")
#
# x_basis_S1_A1 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 1, W, A = 1)))
# x_basis_S0_A1 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 0, W, A = 1)))
# x_basis_S1_A0 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 1, W, A = 0)))
# x_basis_S0_A0 <- make_counter_design_matrix(fit_Q_bar$basis_list, as.matrix(data.table(S = 0, W, A = 0)))
#
# Q_bar_S1_A1_pred <- x_basis_S1_A1 %*% fit_Q_bar$beta
# Q_bar_S0_A1_pred <- x_basis_S0_A1 %*% fit_Q_bar$beta
# Q_bar_S1_A0_pred <- x_basis_S1_A0 %*% fit_Q_bar$beta
# Q_bar_S0_A0_pred <- x_basis_S0_A0 %*% fit_Q_bar$beta
#
# tau_A1_pred <- Q_bar_S1_A1_pred - Q_bar_S0_A1_pred
# tau_A0_pred <- Q_bar_S1_A0_pred - Q_bar_S0_A0_pred
#
# # target parameter point estimate
# psi_pound_pred <- mean((1 - Pi_A0_pred) * tau_A0_pred - (1 - Pi_A1_pred) * tau_A1_pred)
# psi_tilde_pred <- mean(Q_bar_S1_A1_pred * Pi_A1_pred - Q_bar_S1_A0_pred * Pi_A0_pred +
#                        Q_bar_S0_A1_pred * (1 - Pi_A1_pred) - Q_bar_S0_A0_pred * (1 - Pi_A0_pred))
# psi_pred <- psi_tilde_pred - psi_pound_pred
# print(psi_pred)
#
# # TODO: inference
