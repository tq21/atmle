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
                       lambda_method = "glm",
                       v_folds = 5,
                       g_bounds = c(0.01, 0.99),
                       G_bar_bounds = c(0.01, 0.99)) {

  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }

  # define Delta_G=I(C>=t0)
  data[, Delta_G := as.numeric(..Delta == 1 | (..Delta == 0 & t0 < ..T_tilde))]

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
               g_bounds = g_bounds)

  # estimate nuisance: \bar{G}(t0|W,A=a)=P(I(C>=t0)|W,A=a)
  G_bar <- learn_G_bar(W = W,
                       A = A,
                       Delta_G = Delta_G,
                       g = g,
                       method = G_bar_method,
                       folds = folds,
                       G_bar_bounds = G_bar_bounds)

  # estimate nuisance: \Tilde{\theta}(W)=E(Y|W)
  theta_tilde <- learn_theta_tilde(data = data,
                                   W = W,
                                   Y = Y,
                                   delta = data[["Delta_G"]],
                                   method = theta_method,
                                   family = family,
                                   theta_bounds = theta_bounds,
                                   folds = folds)

  # learn a working model for difference in conditional survival functions
  stablize_weights <- g*(1-g)*G_bar$integrate_A
  cate_surv <- learn_T(data = data,
                       W = W,
                       A = A,
                       Y = Y,
                       delta = data[["Delta_G"]],
                       g = g,
                       theta_tilde = theta_tilde,
                       weights = stablize_weights,
                       method = cate_surv_working_model,
                       folds = folds,
                       enumerate_basis_args = enumerate_basis_args,
                       fit_hal_args = fit_hal_args)

  # make long format of data
  data_long <- copy(data)
  n <- data_long[, .N]
  data_long <- data_long[rep(1:.N, each = tau)]
  data_long[, `:=` (id = rep(seq(n), each = tau),
                    t = rep(seq(tau), n))]
  data_long[, T_tilde_t := as.numeric(..T_tilde == t)]

  # make long format of data (failure times only)
  data_fail <- data[Delta_G == 1]
  n_fail <- data_fail[, .N]
  data_fail_long <- copy(data_fail)
  data_fail_long <- data_fail_long[rep(1:.N, each = tau)]
  data_fail_long[, `:=` (id = rep(seq(n_fail), each = tau),
                         t = rep(seq(tau), n_fail))]
  data_fail_long[, T_tilde_t := as.numeric(..T_tilde == t)]

  # estimate nuisance: lambda(t|W,A)=P(T=t|T>=t,W,A)
  lambda <- learn_lambda(data_fail_long = data_fail_long,
                         data_long = data_long,
                         W = W,
                         A = A,
                         T_tilde = T_tilde,
                         Delta = Delta,
                         tau = tau,
                         method = lambda_method,
                         folds = folds)











  # compute survival probabilities from hazard
  data_long[, `:=` (lambda = lambda$pred,
                    lambda_A1 = lambda$A1,
                    lambda_A0 = lambda$A0)]
  data_long[, `:=` (surv_A1 = cumprod(1 - lambda_A1),
                    surv_A0 = cumprod(1 - lambda_A0)), by = id]

  # compute CATE=S(t0|A=1,W)-S(t0|A=0,W)
  Q_bar_r <- lambda[t == t0, surv_A1] - lambda[t == t0, surv_A0]

  # learn working model for conditional effect of A
  Q_bar_r_working_model <- learn_Q_bar_r_working(W = W,
                                                 Q_bar_r = Q_bar_r,
                                                 G_bar = G_bar,
                                                 weights = Delta_G / G_bar$integrate_A)

  # tmle targeting of lambda (this requires a working model for Q_bar_r)
  lambda <- lambda_tmle(A = A,
                        T_tilde = T_tilde,
                        t0 = t0,
                        g = g,
                        G_bar = G_bar,
                        lambda = lambda,
                        Q_bar_r_working_model = Q_bar_r_working_model)

  # compute efficient influence curve
  psi_tilde_eic <- get_eic_psi_tilde_surv(g = g,
                                          G_bar = G_bar,
                                          Q_bar_r_working_model = Q_bar_r_working_model,
                                          Q_bar_r = Q_bar_r,
                                          lambda = lambda,
                                          t_seq = sort(unique(T_tilde)))

  psi_tilde_est <- mean(Q_bar_r)
  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
  psi_tilde_lower <- psi_tilde_est - 1.96 * psi_tilde_se
  psi_tilde_upper <- psi_tilde_est + 1.96 * psi_tilde_se

  return(list(psi_tilde = psi_tilde_est,
              psi_tilde_lower = psi_tilde_lower,
              psi_tilde_upper = psi_tilde_upper))
}
