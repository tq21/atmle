# data structure: O=(S,W,A,T_tilde=min(T,C),Delta=I(T<=C))
atmle_surv <- function(data,
                       S,
                       W,
                       A,
                       T_tilde,
                       Delta,
                       t0,
                       g_rct,
                       controls_only,
                       g_method = "glmnet",
                       G_bar_method = "glmnet",
                       theta_method = "glmnet",
                       lambda_method = "glm",
                       cate_surv_working_model = "glmnet",
                       v_folds = 5,
                       g_bounds = c(0.01, 0.99),
                       G_bar_bounds = c(0.01, 0.99),
                       theta_bounds = c(-Inf, Inf)) {

  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  tau <- max(data[[T_tilde]])

  # define Delta_G=I(C>=t0), Y=I(T>=t0)
  set(data, j = "Delta_G", value = as.numeric(data[[Delta]] == 1 | (data[[Delta]] == 0 & data[[T_tilde]] < t0)))
  set(data, j = "Y", value = as.numeric(data[[Delta]] == 1 & data[[T_tilde]] >= t0))

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(data[, ..W]) == 1) {
    data[, W_dummy := 1]
    W <- c(W, "W_dummy")
  }

  folds <- make_folds(n = n, V = v_folds)

  # psi tilde ------------------------------------------------------------------
  # estimate nuisance: g(1|W)=P(A=1|W)
  g <- learn_g(data = data,
               W = W,
               A = A,
               method = g_method,
               folds = folds,
               g_bounds = g_bounds) # GOOD!

  # estimate nuisance: \bar{G}(t0|W,A=a)=P(I(C>=t0)|W,A=a)
  G_bar <- learn_G_bar(data = data,
                       W = W,
                       A = A,
                       g = g,
                       method = G_bar_method,
                       folds = folds,
                       G_bar_bounds = G_bar_bounds) # GOOD!

  # estimate nuisance: \Tilde{\theta}(W)=E(Y|W)
  theta_tilde <- learn_theta_tilde(data = data,
                                   W = W,
                                   Y = "Y",
                                   delta = data[["Delta_G"]],
                                   method = theta_method,
                                   family = "binomial",
                                   theta_bounds = theta_bounds,
                                   folds = folds) # GOOD!

  # learn a working model for difference in conditional survival functions
  stablize_weights <- g*(1-g)*G_bar$integrate_A
  cate_surv <- learn_T(data = data,
                       W = W,
                       A = A,
                       Y = "Y",
                       delta = data[["Delta_G"]],
                       g = g,
                       theta_tilde = theta_tilde,
                       weights = stablize_weights,
                       method = cate_surv_working_model,
                       folds = folds,
                       enumerate_basis_args = enumerate_basis_args,
                       fit_hal_args = fit_hal_args) # GOOD!

  # make long format of data
  data_long <- copy(data)
  n <- data_long[, .N]
  data_long <- data_long[rep(1:.N, each = tau)]
  data_long[, `:=` (id = rep(seq(n), each = tau),
                    t = rep(seq(tau), n))]
  data_long[, T_t := as.numeric((T_tilde) == t & (Delta) == 1)]

  # # make long format of data (failure times only)
  # data_fail <- data[Delta_G == 1]
  # n_fail <- data_fail[, .N]
  # data_fail_long <- copy(data_fail)
  # data_fail_long <- data_fail_long[rep(1:.N, each = tau)]
  # data_fail_long[, `:=` (id = rep(seq(n_fail), each = tau),
  #                        t = rep(seq(tau), n_fail))]
  # set(data_fail_long, j = "T_tilde_t", value = as.numeric(data_fail_long[[T_tilde]] == data_fail_long[["t"]]))

  # estimate nuisance: lambda(t|W,A)=P(T=t|T>=t,W,A)
  lambda <- learn_lambda(data_long = data_long,
                         W = W,
                         A = A,
                         T_tilde = T_tilde,
                         method = lambda_method,
                         folds = folds) # GOOD!

  # compute survival probabilities from hazard
  data_long[, `:=` (lambda = lambda$pred,
                    lambda_A1 = lambda$A1,
                    lambda_A0 = lambda$A0)]
  data_long[, `:=` (surv_A1 = cumprod(1 - lambda_A1),
                    surv_A0 = cumprod(1 - lambda_A0)), by = id]

  # TMLE targeting of lambda
  data_long <- lambda_tmle(data_long = data_long,
                           n = data[, .N],
                           A = A,
                           T_tilde = T_tilde,
                           t0 = t0,
                           g = g,
                           G_bar = G_bar,
                           lambda = lambda,
                           stablize_weights = stablize_weights,
                           cate_surv = cate_surv) # GOOD!

  return(list(psi_tilde_est = mean(cate_surv$pred)))

  # Re-learn beta under the survival mapped from the targeted lambda
  data_long[t == t0-1, `:=` (targeted_diff = surv_A1 - surv_A0)]
  targeted_diff <- data_long[t == t0-1, targeted_diff]
  targeted_working_model <- glm(targeted_diff ~ -1+.,
                                data = as.data.frame(cate_surv$x_basis),
                                family = "gaussian") # TODO: can bound this using logistic
  cate_surv$coefs <- as.numeric(coef(targeted_working_model))
  cate_surv$coefs[is.na(cate_surv$coefs)] <- 0
  cate_surv$pred <- as.numeric(cate_surv$x_basis %*% cate_surv$coefs)

  # unique_t <- sort(unique(data[[T_tilde]]))

  psi_tilde_est <- mean(cate_surv$pred)

  # psi_tilde_est <- mean(data_long[t == t0-1, surv_A1] - data_long[t == t0-1, surv_A0])
  return(list(psi_tilde_est = psi_tilde_est))

  # compute efficient influence curve
  psi_tilde_eic <- get_eic_psi_tilde_surv(data = data,
                                          data_long = data_long,
                                          g = g,
                                          stablize_weight = stablize_weights,
                                          cate_surv = cate_surv,
                                          unique_t = unique_t,
                                          Y = "Y",
                                          n = data[, .N])

  # psi_tilde_est <- mean(Q_bar_r)
  # psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
  # psi_tilde_lower <- psi_tilde_est - 1.96 * psi_tilde_se
  # psi_tilde_upper <- psi_tilde_est + 1.96 * psi_tilde_se
#
  # return(list(psi_tilde = psi_tilde_est,
  #             psi_tilde_lower = psi_tilde_lower,
  #             psi_tilde_upper = psi_tilde_upper))
}
