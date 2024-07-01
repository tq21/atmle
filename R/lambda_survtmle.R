lambda_survtmle <- function(data_long,
                            n,
                            A,
                            T_tilde,
                            g,
                            t0,
                            lambda) {

  unique_t <- sort(unique(data_long[[T_tilde]]))

  # clever covariate
  data_long[, surv_t0_A1 := surv_A1[t == t0], by = id]
  data_long[, surv_t0_A0 := surv_A0[t == t0], by = id]
  data_long[, H_A1 := -((A)/(g*surv_c_A1_lag)*(surv_t0_A1/surv_A1))*(t<=t0)]
  data_long[, H_A0 := (-(1-(A))/((1-g)*surv_c_A0_lag)*(surv_t0_A0/surv_A0))*(t<=t0)]
  data_long[, clever_cov_A1 := as.numeric((T_tilde) >= t)*H_A1]
  data_long[, clever_cov_A0 := as.numeric((T_tilde) >= t)*H_A0]
  data_long[, dN_t := as.numeric((T_tilde) == t & (Delta) == 1)]

  # TMLE: local least favorable submodel update, from tau (last time point)
  unique_t_rev <- seq(t0, 1)
  for (j in 1:length(unique_t_rev)) {
    cur_t <- unique_t_rev[j]
    data_sub <- data_long[t == cur_t]
    epsilons <- coef(glm(dN_t ~ -1 + offset(qlogis(lambda)) + clever_cov_A1 + clever_cov_A0,
                         data = data_sub, family = "quasibinomial"))
    epsilons[is.na(epsilons)] <- 0

    # tmle update
    data_long[t == cur_t, lambda_A1 := plogis(qlogis(data_sub$lambda_A1) + epsilons[1] * data_sub$clever_cov_A1)]
    data_long[t == cur_t, lambda_A0 := plogis(qlogis(data_sub$lambda_A0) + epsilons[2] * data_sub$clever_cov_A0)]
    data_long[t == cur_t, lambda := plogis(qlogis(data_sub$lambda) + epsilons[1] * data_sub$clever_cov_A1 + epsilons[2] * data_sub$clever_cov_A0)]

    # recompute relevant parts using updated lambda
    data_long[, `:=` (surv_A1 = cumprod(1 - lambda_A1),
                      surv_A0 = cumprod(1 - lambda_A0)), by = id]
    data_long[, surv_t0_A1 := surv_A1[t == t0], by = id]
    data_long[, surv_t0_A0 := surv_A0[t == t0], by = id]
    data_long[, H_A1 := ((A)/(g*surv_c_A1_lag)*(surv_t0_A1/surv_A1))*(t<=t0)]
    data_long[, H_A0 := (-(1-(A))/((1-g)*surv_c_A0_lag)*(surv_t0_A0/surv_A0))*(t<=t0)]
    data_long[, clever_cov_A1 := as.numeric((T_tilde) >= t)*H_A1]
    data_long[, clever_cov_A0 := as.numeric((T_tilde) >= t)*H_A0]
  }

  return(data_long)
}
