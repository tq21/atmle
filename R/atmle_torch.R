library(torch)
atmle_torch <- function(data,
                        S,
                        W,
                        A,
                        Y,
                        controls_only,
                        family,
                        eic_method,
                        target_gwt = TRUE,
                        lr = 1e-2,
                        max_iter = 5000,
                        patience = 10,
                        tolerance = 1e-6,
                        theta_method = "glm",
                        theta_tilde_method = "glm",
                        g_method = "glm",
                        g_bounds = c(0.01, 0.99),
                        Pi_bounds = c(0.01, 0.99),
                        theta_bounds = c(-Inf, Inf),
                        enumerate_basis_args = list(max_degree = 2,
                                                    smoothness_orders = 1),
                        v_folds = 10,
                        parallel = FALSE,
                        cross_fit_nuisance = TRUE,
                        verbose = FALSE,
                        device = "cpu",
                        browse = FALSE) {

  if (browse) browser()

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # extract variables ----------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  delta <- as.integer(!is.na(Y))
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
  if (family == "gaussian") {
    if (sum(delta) < n) {
      cv_strata <- paste0(S, "-", A, "-", delta)
    } else {
      cv_strata <- paste0(S, "-", A)
    }
  } else if (family == "binomial") {
    if (sum(delta) < n) {
      cv_strata <- paste0(S, "-", A, "-", Y, "-", delta)
    } else {
      cv_strata <- paste0(S, "-", A, "-", Y)
    }
  }
  suppressWarnings({
    folds <- make_folds(
      n = n, V = v_folds,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # learn nuisance parameters
  theta_tilde <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = delta,
    method = theta_tilde_method,
    folds = folds,
    family = family,
    theta_bounds = theta_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  theta <- learn_theta(
    W = W,
    A = A,
    Y = Y,
    delta = delta,
    controls_only = controls_only,
    method = theta_method,
    folds = folds,
    family = family,
    theta_bounds = theta_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  g <- learn_g(
    S = S,
    W = W,
    A = A,
    method = g_method,
    controls_only = controls_only,
    v_folds = v_folds,
    g_bounds = g_bounds
  )

  Pi <- learn_Pi(
    g = g,
    A = A,
    Pi_bounds = Pi_bounds
  )

  if (sum(delta) < n) {
    # outcome has missing
    g_delta <- learn_g_delta(
      S = S,
      W = W,
      A = A,
      delta = delta,
      method = g_delta_method,
      folds = folds,
      g_bounds = g_bounds
    )
    if (verbose) cat("Done!\n")

    g_delta_tilde <- learn_g_delta_tilde(
      W = W,
      A = A,
      delta = delta,
      method = g_delta_method,
      folds = folds,
      g_bounds = g_bounds
    )
  } else {
    # no censoring
    g_delta <- g_delta_tilde <- list(
      pred = rep(1, length(A)),
      A0 = rep(1, length(A)),
      A1 = rep(1, length(A))
    )
  }

  # censoring weights
  weights <- delta / g_delta$pred
  weights_tilde <- delta / g_delta_tilde$pred

  # learn CATE
  tau_A_seq <- learn_tau_A_seq(
    S = S,
    W = W,
    A = A,
    Y = Y,
    theta = theta_tilde,
    Q_bar = theta,
    g1W = g$pred,
    Pi1WA = Pi$pred,
    weights = weights,
    enumerate_basis_args = enumerate_basis_args,
    lr = lr,
    max_iter = max_iter,
    verbose = verbose,
    device = device,
    tolerance = tolerance,
    patience = patience,
    parallel = parallel
  )

  # learn CARE
  tau_S_seq <- learn_tau_S_seq(
    S = S,
    W = W,
    A = A,
    Y = Y,
    theta = theta_tilde,
    Q_bar = theta,
    g1W = g$pred,
    Pi1WA = Pi$pred,
    weights = weights,
    enumerate_basis_args = enumerate_basis_args,
    lr = lr,
    max_iter = max_iter,
    verbose = verbose,
    device = device,
    tolerance = tolerance,
    patience = patience,
    parallel = parallel
  )

  # TMLE to target Pi
  min_len <- min(length(tau_A_seq$beta_list), length(tau_S_seq$beta_list))
  res_seq <- map(seq(min_len), function(j) {
    tau_A_cur <- tau_A_seq$beta_list[[j]]
    tau_S_cur <- tau_S_seq$beta_list[[j]]

    tau_S_tmp <- vector("list", length = 6)
    tau_S_tmp$x_basis <- as.matrix(cbind(1, tau_S_seq$hal_design[, tau_S_cur$idx, drop = FALSE]))
    tau_S_tmp$x_basis_A1 <- as.matrix(cbind(1, tau_S_seq$hal_design_A1[, tau_S_cur$idx, drop = FALSE]))
    tau_S_tmp$x_basis_A0 <- as.matrix(cbind(1, tau_S_seq$hal_design_A0[, tau_S_cur$idx, drop = FALSE]))
    tau_S_tmp$A1 <- as.vector(tau_S_tmp$x_basis_A1 %*% as.matrix(tau_S_cur$beta))
    tau_S_tmp$A0 <- as.vector(tau_S_tmp$x_basis_A0 %*% as.matrix(tau_S_cur$beta))
    tau_S_tmp$pred <- as.vector(tau_S_tmp$x_basis %*% as.matrix(tau_S_cur$beta))

    tau_A_tmp <- vector("list", length = 2)
    tau_A_tmp$x_basis <- as.matrix(cbind(1, tau_A_seq$hal_design[, tau_A_cur$idx, drop = FALSE]))
    tau_A_tmp$pred <- as.vector(tau_A_tmp$x_basis %*% as.matrix(tau_A_cur$beta))

    # target Pi
    Pi_cur <- Pi_tmle(
      S = S,
      W = W,
      A = A,
      g = g$pred,
      tau = tau_S_tmp,
      Pi = Pi,
      controls_only = controls_only,
      target_gwt = target_gwt,
      Pi_bounds = Pi_bounds
    )

    # point estimates
    psi_tilde_est <- mean(tau_A_tmp$pred)
    if (controls_only) {
      psi_pound_est <- mean((1 - Pi_cur$A0) * tau_S_tmp$A0)
    } else {
      psi_pound_est <- mean((1 - Pi_cur$A0) * tau_S_tmp$A0 - (1 - Pi_cur$A1) * tau_S_tmp$A1)
    }

    # inference
    psi_tilde_eic <- get_eic_psi_tilde(
      psi_tilde = tau_A_tmp,
      g = g$pred,
      theta = theta_tilde,
      Y = Y,
      A = A,
      n = n,
      weights = weights_tilde,
      eic_method = eic_method
    )
    psi_pound_eic <- get_eic_psi_pound(
      Pi = Pi_cur,
      tau = tau_S_tmp,
      g = g$pred,
      theta = theta,
      psi_pound_est = psi_pound_est,
      S = S,
      A = A,
      Y = Y,
      n = n,
      controls_only = controls_only,
      weights = weights,
      eic_method = eic_method
    )
    psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE) / n)
    psi_tilde_lower <- psi_tilde_est - 1.96 * psi_tilde_se
    psi_tilde_upper <- psi_tilde_est + 1.96 * psi_tilde_se
    psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE) / n)
    psi_pound_lower <- psi_pound_est - 1.96 * psi_pound_se
    psi_pound_upper <- psi_pound_est + 1.96 * psi_pound_se
    est <- psi_tilde_est - psi_pound_est
    eic <- psi_tilde_eic - psi_pound_eic
    se <- sqrt(var(eic, na.rm = TRUE) / n)
    lower <- est - 1.96 * se
    upper <- est + 1.96 * se
    sn <- sqrt(var(eic, na.rm = TRUE))/(sqrt(n) * log(n))

    return(list(psi = est,
                lower = lower,
                upper = upper,
                PnEIC = mean(eic),
                sn = sn,
                tau_S_beta = tau_S_cur$beta,
                tau_A_beta = tau_A_cur$beta))
  })

  return(res_seq)
}
