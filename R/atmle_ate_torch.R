library(torch)
atmle_ate_torch <- function(data,
                            W,
                            A,
                            Y,
                            family,
                            eic_method,
                            lr = 1e-2,
                            max_iter = 5000,
                            patience = 10,
                            tolerance = 1e-6,
                            theta_method = "glm",
                            g_method = "glm",
                            g_bounds = c(0.01, 0.99),
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
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  n <- nrow(data)

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(W) == 1) {
    W <- cbind(W, 1)
  }

  # cross fitting schemes
  if (family == "gaussian") {
    folds <- make_folds(n = n, V = v_folds)
  } else if (family == "binomial") {
    folds <- make_folds(n = n, V = v_folds, strata_ids = as.integer(factor(Y)))
  }

  # learn nuisance parameters
  theta <- learn_theta_tilde(
    W = W,
    Y = Y,
    delta = rep(1, n),
    method = theta_method,
    folds = folds,
    family = family,
    theta_bounds = theta_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  g <- learn_g_tmp(
    W = W,
    A = A,
    method = g_method,
    folds = folds,
    g_bounds = g_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  # learn CATE
  tau_A_seq <- learn_tau_A_seq(
    S = rep(NA, n),
    W = W,
    A = A,
    Y = Y,
    theta = theta,
    Q_bar = rep(NA, n),
    g1W = g,
    Pi1WA = rep(NA, n),
    weights = rep(1, n),
    enumerate_basis_args = enumerate_basis_args,
    lr = lr,
    max_iter = max_iter,
    verbose = verbose,
    device = device,
    tolerance = tolerance,
    patience = patience,
    parallel = parallel
  )

  # estimator sequence
  res_seq <- map(seq(length(tau_A_seq$beta_list)), function(j) {
    tau_A_cur <- tau_A_seq$beta_list[[j]]
    tau_A_tmp <- vector("list", length = 2)
    tau_A_tmp$x_basis <- as.matrix(cbind(1, tau_A_seq$hal_design[, tau_A_cur$idx, drop = FALSE]))
    tau_A_tmp$pred <- as.vector(tau_A_tmp$x_basis %*% as.matrix(tau_A_cur$beta))

    # point estimates and inference
    psi <- mean(tau_A_tmp$pred)
    eic <- get_eic_psi_tilde(
      psi_tilde = tau_A_tmp,
      g = g,
      theta = theta,
      Y = Y,
      A = A,
      n = n,
      weights = rep(1, n),
      eic_method = eic_method
    )
    se <- sqrt(var(eic, na.rm = TRUE) / n)
    lower <- psi - 1.96 * se
    upper <- psi + 1.96 * se
    sn <- sqrt(var(eic, na.rm = TRUE))/(sqrt(n) * log(n))

    return(list(psi = psi,
                lower = lower,
                upper = upper,
                EIC = eic,
                PnEIC = mean(eic),
                sn = sn,
                L1_norm = sum(abs(tau_A_cur$beta))))
  })

  res_seq <- res_seq[sort(map_vec(res_seq, function(.x) .x$L1_norm), index.return = TRUE)$ix]

  return(res_seq)
}
