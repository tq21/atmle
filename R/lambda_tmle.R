lambda_tmle <- function(A,
                        T_tilde,
                        t0,
                        g,
                        G_bar,
                        lambda,
                        Q_bar_r_working_model) {

  unique_T_tilde <- sort(unique(T_tilde))
  t0_idx <- which(t0 == unique_T_tilde)
  t_idx_grid <- seq(t0_idx)

  # clever covariate
  lambda[, surv_t0_A1 := surv_A1[t == t0], by = id]
  lambda[, surv_t0_A0 := surv_A0[t == t0], by = id]

  IM <- -t(Q_bar_r_working_model$phi) %*% diag(g*(1-g)) %*% Q_bar_r_working_model$phi / length(A)
  comp_1 <- as.numeric(Q_bar_r_working_model$phi %*% solve(IM) %*% colMeans(Q_bar_r_working_model$phi) * (g * (1 - g) * G_bar$integrate_A))

  id_to_nuisance <- data.table(id = seq(length(A)),
                               A = A,
                               T_tilde = T_tilde,
                               g = g,
                               G_bar_A1 = G_bar$A1,
                               G_bar_A0 = G_bar$A0,
                               comp_1 = comp_1)
  lambda <- merge(lambda, id_to_nuisance, by = "id", all.x = TRUE) # TODO: this step sorted the rows in a different way, avoid that
  lambda[, `:=` (H_A1 = -A/(g*G_bar_A1)*(surv_t0_A1/surv_A1),
                 H_A0 = -(1-A)/((1-g)*G_bar_A0)*(surv_t0_A0/surv_A0),
                 T_tilde_geq_t = as.numeric(T_tilde >= t))]
  lambda[, `:=` (clever_cov_A1 = comp_1 * T_tilde_geq_t * H_A1,
                 clever_cov_A0 = comp_1 * T_tilde_geq_t * (-H_A0),
                 dN_t = as.numeric(T_tilde == t))] # TODO: check definition of dN_t

  # TODO: TMLE: local least favorable submodel update, from tau (last time point)
  all_t_rev <- rev(unique_T_tilde)

  for (j in 1:length(all_t_rev)) {
    cur_t <- all_t_rev[j]
    sub_lambda <- lambda[t == cur_t]
    sub_lambda[, lambda_pred := A * lambda_A1 + (1 - A) * lambda_A0]
    epsilons <- coef(glm(dN_t ~ -1 + offset(qlogis(lambda_pred)) + clever_cov_A1 + clever_cov_A0,
                         data = sub_lambda, family = "quasibinomial"))
    epsilons[is.na(epsilons)] <- 0

    # tmle update
    lambda[t == cur_t, lambda_A1 := plogis(qlogis(sub_lambda$lambda_A1) + epsilons[1] * sub_lambda$clever_cov_A1)]
    lambda[t == cur_t, lambda_A0 := plogis(qlogis(sub_lambda$lambda_A0) + epsilons[2] * sub_lambda$clever_cov_A0)]
    lambda[t == cur_t, lambda_pred := plogis(
      qlogis(sub_lambda$lambda_pred) + epsilons[1] * clever_cov_A1 + epsilons[2] * clever_cov_A0)]

    # TODO: make sure keys are id and t
  }

  return(lambda)
}
