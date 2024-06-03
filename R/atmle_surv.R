# data structure: O=(S,W,A,T_tilde=min(T,C),Delta=I(T<=C))
atmle_surv <- function(data,
                       S_node,
                       W_node,
                       A_node,
                       T_tilde_node,
                       Delta_node,
                       t0,
                       g_rct,
                       controls_only,
                       g_method = "glmnet",
                       G_bar_method = "glmnet",
                       lambda_method = "glmnet",
                       v_folds = 5,
                       g_bounds = c(0.01, 0.99),
                       G_bar_bounds = c(0.01, 0.99)) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # define nodes ---------------------------------------------------------------
  S <- data[, S_node]
  W <- data[, W_node, drop = FALSE]
  A <- data[, A_node]
  T_tilde <- data[, T_tilde_node]
  Delta <- data[, Delta_node]
  n <- nrow(data)

  # validate controls_only argument
  if (controls_only & 1 %in% A[S == 0]) {
    stop("The 'controls_only' argument is set to TRUE, but there are treated units in the external data.")
  } else if (!controls_only & sum(A[S == 0]) == 0) {
    stop("The 'controls_only' argument is set to FALSE, but there are only control units in the external data.")
  }

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(W) == 1) {
    W <- cbind(W, 1)
  }

  # cross fitting schemes
  # cv_strata <- paste0(S, "-", A)
  # suppressWarnings({
  #   folds <- make_folds(
  #     n = n, V = v_folds,
  #     strata_ids = as.integer(factor(cv_strata))
  #   )
  # })

  folds <- make_folds(n = n, V = v_folds)

  # psi tilde ------------------------------------------------------------------
  # estimate nuisance: g(1|W)=P(A=1|W)
  g <- learn_g(W = W,
               A = A,
               method = g_method,
               folds = folds,
               g_bounds = g_bounds)

  # estimate nuisance: \bar{G}(t0|W,A=a)=P(I(C>=t0)|W,A=a)
  cens_after_t0 <- as.numeric(Delta == 0 & T_tilde >= t0)
  G_bar <- learn_G_bar(W = W,
                       A = A,
                       cens_after_t0 = cens_after_t0,
                       g = g,
                       method = G_bar_method,
                       folds = folds,
                       G_bar_bounds = G_bar_bounds)

  # estimate nuisance: lambda(t|W,A)=P(T=t|T>=t,W,A)
  lambda <- learn_lambda(W = W,
                         A = A,
                         T_tilde = T_tilde,
                         method = lambda_method,
                         folds = folds)

  # compute survival probabilities from hazard
  lambda[, `:=` (surv_A1 = cumprod(1 - lambda_A1),
                 surv_A0 = cumprod(1 - lambda_A0)), by = id]

  # compute CATE=S(t0|A=1,W)-S(t0|A=0,W)
  Q_bar_r <- lambda[t == t0, surv_A1] - lambda[t == t0, surv_A0]

  # learn working model for conditional effect of A
  Q_bar_r_working_model <- learn_Q_bar_r_working(W = W,
                                                 Q_bar_r = Q_bar_r,
                                                 G_bar = G_bar)

  # tmle targeting of lambda (this requires a working model for Q_bar_r)
  lambda <- lambda_tmle(A = A,
                        T_tilde = T_tilde,
                        t0 = t0,
                        g = g,
                        G_bar = G_bar,
                        lambda = lambda,
                        Q_bar_r_working_model = Q_bar_r_working_model)

  # compute efficient influence curve


}






