lambda_tmle <- function(data_long,
                        n,
                        A,
                        T_tilde,
                        t0,
                        g,
                        G_bar,
                        lambda,
                        stablize_weights,
                        cate_surv) {

  unique_t <- sort(unique(data_long[[T_tilde]]))

  # clever covariate
  data_long[, surv_t0_A1 := surv_A1[t == t0], by = id]
  data_long[, surv_t0_A0 := surv_A0[t == t0], by = id]
  data_long[, H_A1 := ((A)/(g*surv_c_A1)*(surv_t0_A1/surv_A1))*(t<=t0)]
  tmp <- (data_long[[A]]/(g*data_long$surv_c_A1)*(data_long$surv_t0_A1/data_long$surv_A1))*(data_long$t<=t0)
  if (any(tmp != data_long$H_A1)) stop("H_A1 is not computed correctly")

  data_long[, H_A0 := (-(1-(A))/((1-g)*surv_c_A0)*(surv_t0_A0/surv_A0))*(t<=t0)]
  tmp <- (-(1-(data_long[[A]]))/((1-g)*data_long$surv_c_A0)*(data_long$surv_t0_A0/data_long$surv_A0))*(data_long$t<=t0)
  if (any(tmp != data_long$H_A0)) stop("H_A0 is not computed correctly")

  IM <- -t(cate_surv$x_basis) %*% diag(g*(1-g)) %*% cate_surv$x_basis / n
  tmp <- as.numeric(cate_surv$x_basis %*% solve(IM) %*% colMeans(cate_surv$x_basis) * stablize_weights)
  data_long[, clever_cov_A1 := tmp*as.numeric((T_tilde) >= t)*H_A1]
  data_long[, clever_cov_A0 := tmp*as.numeric((T_tilde) >= t)*H_A0]
  data_long[, dN_t := as.numeric((T_tilde) == t & Delta_G == 1)]

  # TMLE: local least favorable submodel update, from tau (last time point)
  unique_t_rev <- rev(unique_t)
  for (j in 1:length(unique_t_rev)) {
    cur_t <- unique_t_rev[j]
    data_sub <- data_long[t == cur_t]
    epsilons <- coef(glm(dN_t ~ -1 + offset(qlogis(lambda)) + clever_cov_A1 + clever_cov_A0,
                         data = data_sub, family = "quasibinomial"))
    epsilons[is.na(epsilons)] <- 0

    # tmle update
    data_long[t == cur_t, lambda_A1 := plogis(qlogis(data_sub$lambda_A1) + epsilons[1] * data_sub$clever_cov_A1)]
    data_long[t == cur_t, lambda_A0 := plogis(qlogis(data_sub$lambda_A0) + epsilons[2] * data_sub$clever_cov_A0)]
    data_long[t == cur_t, lambda := plogis(qlogis(data_sub$lambda) + epsilons[1] * clever_cov_A1 + epsilons[2] * clever_cov_A0)]
  }

  # recompute relevant parts using updated lambda
  data_long[, `:=` (surv_A1 = shift(cumprod(1 - lambda_A1), fill = 1),
                    surv_A0 = shift(cumprod(1 - lambda_A0), fill = 1)), by = id]
  data_long[, surv_t0_A1 := surv_A1[t == t0], by = id]
  data_long[, surv_t0_A0 := surv_A0[t == t0], by = id]
  data_long[, H_A1 := ((A)/(g*G_bar$A1)*(surv_t0_A1/surv_A1))*(t<=t0)]
  data_long[, H_A0 := (-(1-(A))/((1-g)*G_bar$A0)*(surv_t0_A0/surv_A0))*(t<=t0)]
  data_long[, clever_cov_A1 := tmp*as.numeric((T_tilde) >= t)*H_A1]
  data_long[, clever_cov_A0 := tmp*as.numeric((T_tilde) >= t)*H_A0]

  return(data_long)
}
