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
                       theta_bounds = c(0.01, 0.99)) {

  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  tau <- max(data[[T_tilde]])

  set(data, j = "Delta_t0", value = as.numeric(data[[Delta]] == 1 | data[[T_tilde]] > t0))
  set(data, j = "Y", value = as.numeric(data[[T_tilde]] > t0))
  set(data, j = "U", value = pmin(data[[T_tilde]], t0))

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(data[, ..W]) == 1) {
    data[, W_dummy := 1]
    W <- c(W, "W_dummy")
  }

  folds <- make_folds(n = n, V = v_folds)

  # make long format of data
  data_long <- copy(data)
  n <- data_long[, .N]
  data_long <- data_long[rep(1:.N, each = tau)]
  data_long[, `:=` (id = rep(seq(n), each = tau),
                    t = rep(seq(tau), n))]
  data_long[, T_t := as.numeric((T_tilde) == t & (Delta) == 1)]
  data_long[, C_t := as.numeric((T_tilde) == t & (Delta) == 0)]
  data_long[, C_tilde := tau][(Delta) == 0, C_tilde := (T_tilde)]

  # psi tilde ------------------------------------------------------------------
  # estimate nuisance: g(1|W)=P(A=1|W)
  g <- learn_g(data = data,
               W = W,
               A = A,
               method = g_method,
               folds = folds,
               g_bounds = g_bounds) # GOOD! 6/23/2024

  # estimate nuisance: \bar{G}(t|W,A)=P(C>t|W,A)
  lambda_c <- learn_hazard(data_long = data_long,
                           X = c(W, A),
                           T_tilde = "C_tilde",
                           event = "C_t",
                           method = lambda_method,
                           folds = folds,
                           counter_var = A) # GOOD! 6/23/2024

  # compute survival function of censoring from hazard estimates
  data_long[, `:=` (lambda_c = lambda_c$pred,
                    lambda_c_A1 = lambda_c$A1,
                    lambda_c_A0 = lambda_c$A0)]
  data_long[, `:=` (surv_c_A1 = cumprod(1 - lambda_c_A1),
                    surv_c_A0 = cumprod(1 - lambda_c_A0)), by = id]

  # compute nuisance: \bar{G}(t0|W,A)=P(C>t0|W,A)
  G_bar <- list(A1 = data_long[t == U, ]$surv_c_A1,
                A0 = data_long[t == U, ]$surv_c_A0)
  G_bar$pred <- data[[A]] * G_bar$A1 + (1 - data[[A]]) * G_bar$A0
  G_bar$integrate_A <- as.numeric(G_bar$A1)*g+as.numeric(G_bar$A0)*(1-g)

  # estimate nuisance: \Tilde{\theta}(W)=P(T>t0|W)
  theta_tilde <- learn_hazard(data_long = data_long,
                              X = W,
                              T_tilde = T_tilde,
                              event = "T_t",
                              method = lambda_method,
                              folds = folds,
                              counter_var = NULL) # GOOD! 6/23/2024

  # learn a working model for difference in conditional survival functions
  stablize_weights <- g*(1-g)*G_bar$integrate_A
  cate_surv <- learn_T(data = data,
                       W = W,
                       A = A,
                       Y = "Y",
                       delta = data[["Delta_t0"]],
                       g = g,
                       theta_tilde = theta_tilde,
                       weights = 1/G_bar$pred,
                       method = cate_surv_working_model,
                       folds = folds,
                       enumerate_basis_args = enumerate_basis_args,
                       fit_hal_args = fit_hal_args) # GOOD! 6/23/2024

  psi_tilde_r_learner <- mean(cate_surv$pred) # TODO: TESTING PURPOSES

  # estimate nuisance: lambda(t|W,A)=P(T=t|T>=t,W,A)
  lambda <- learn_hazard(data_long = data_long,
                         X = c(W, A),
                         T_tilde = T_tilde,
                         event = "T_t",
                         method = lambda_method,
                         folds = folds,
                         counter_var = A) # GOOD! 6/23/2024

  # compute survival probabilities from hazard
  data_long[, `:=` (lambda = lambda$pred,
                    lambda_A1 = lambda$A1,
                    lambda_A0 = lambda$A0)]
  data_long[, `:=` (surv_A1 = shift(cumprod(1 - lambda_A1), fill = 1),
                    surv_A0 = shift(cumprod(1 - lambda_A0), fill = 1)), by = id]
  psi_tilde_no_tmle_lambda <- mean(data_long[t == t0, surv_A1]-data_long[t == t0, surv_A0])

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
                           cate_surv = cate_surv)
  psi_tilde_tmle_lambda <- mean(data_long[t == t0, surv_A1]-data_long[t == t0, surv_A0])

  # Re-learn beta under the survival mapped from the targeted lambda
  data_long[t == t0, `:=` (targeted_diff = surv_A1 - surv_A0)]
  targeted_diff <- data_long[t == t0, targeted_diff]
  targeted_working_model <- glm(targeted_diff ~ -1+.,
                                data = as.data.frame(cate_surv$x_basis),
                                family = "gaussian") # TODO: can bound this using logistic
  cate_surv$coefs <- as.numeric(coef(targeted_working_model))
  cate_surv$coefs[is.na(cate_surv$coefs)] <- 0
  cate_surv$pred <- as.numeric(cate_surv$x_basis %*% cate_surv$coefs)
  psi_tilde_r_learner_tmle_beta <- mean(cate_surv$pred)

  # compute efficient influence curve
  unique_t <- sort(unique(data[[T_tilde]]))
  psi_tilde_eic <- get_eic_psi_tilde_surv(data = data,
                                          data_long = data_long,
                                          g = g,
                                          Q_bar_r = targeted_diff,
                                          stablize_weight = stablize_weights,
                                          cate_surv = cate_surv,
                                          unique_t = unique_t,
                                          Y = "Y",
                                          n = data[, .N])

  # psi_tilde_est <- mean(data_long[t == t0-1, surv_A1] - data_long[t == t0-1, surv_A0])

  # psi_tilde_est <- mean(Q_bar_r)
  # psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
  # psi_tilde_lower <- psi_tilde_est - 1.96 * psi_tilde_se
  # psi_tilde_upper <- psi_tilde_est + 1.96 * psi_tilde_se
#
  # return(list(psi_tilde = psi_tilde_est,
  #             psi_tilde_lower = psi_tilde_lower,
  #             psi_tilde_upper = psi_tilde_upper))

  return(list(psi_tilde_r_learner = psi_tilde_r_learner,
              psi_tilde_no_tmle_lambda = psi_tilde_no_tmle_lambda,
              psi_tilde_tmle_lambda = psi_tilde_tmle_lambda,
              psi_tilde_r_learner_tmle_beta = psi_tilde_r_learner_tmle_beta))
}
