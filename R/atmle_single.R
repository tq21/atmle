atmle_single <- function(data,
                         W,
                         A,
                         Y,
                         family,
                         theta_method = "glm",
                         g_method = "glm",
                         tau_A_method = "glmnet",
                         g_bounds = c(0.01, 0.99),
                         theta_bounds = c(-Inf, Inf),
                         enumerate_basis_args = list(max_degree = 2,
                                                     smoothness_orders = 1),
                         v_folds = 10,
                         parallel = TRUE,
                         cross_fit_nuisance = TRUE,
                         verbose = FALSE,
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

  g1W <- learn_g_tmp(
    W = W,
    A = A,
    method = g_method,
    folds = folds,
    g_bounds = g_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  # learn CATE
  tau_A <- learn_tau_A(W = W,
                       A = A,
                       Y = Y,
                       theta = theta,
                       g1W = g1W,
                       delta = rep(1, n),
                       v_folds = v_folds,
                       weights = rep(1, n),
                       enumerate_basis_args = enumerate_basis_args,
                       dx = 1e-4,
                       max_iter = 500)
  tau_A_reg <- tau_A$reg
  tau_A_relax <- tau_A$relax

  # point estimates and inference
  psi_reg <- mean(tau_A_reg$pred)
  psi_relax <- mean(tau_A_relax$pred)
  if (ncol(tau_A_reg$x_basis) == 1) {
    eic_proj_reg <- get_eic_psi_tilde(
      psi_tilde = tau_A_reg,
      g = g1W,
      theta = theta,
      Y = Y,
      A = A,
      n = n,
      weights = rep(1, n),
      eic_method = "svd_pseudo_inv"
    )
    eic_proj_relax <- get_eic_psi_tilde(
      psi_tilde = tau_A_relax,
      g = g1W,
      theta = theta,
      Y = Y,
      A = A,
      n = n,
      weights = rep(1, n),
      eic_method = "svd_pseudo_inv"
    )
  } else {
    eic_proj_reg <- grad_proj(W = W,
                              A = A,
                              Y = Y,
                              theta = theta,
                              g1W = g1W,
                              phi_W = tau_A_reg$x_basis,
                              beta = tau_A_reg$coefs)
    eic_proj_relax <- grad_proj(W = W,
                                A = A,
                                Y = Y,
                                theta = theta,
                                g1W = g1W,
                                phi_W = tau_A_relax$x_basis,
                                beta = tau_A_relax$coefs)
  }

  eic_delta_reg <- get_eic_psi_tilde(
    psi_tilde = tau_A_reg,
    g = g1W,
    theta = theta,
    Y = Y,
    A = A,
    n = n,
    weights = rep(1, n),
    eic_method = "svd_pseudo_inv"
  )
  eic_delta_relax <- get_eic_psi_tilde(
    psi_tilde = tau_A_relax,
    g = g1W,
    theta = theta,
    Y = Y,
    A = A,
    n = n,
    weights = rep(1, n),
    eic_method = "svd_pseudo_inv"
  )

  se_proj_reg <- sqrt(var(eic_proj_reg, na.rm = TRUE) / n)
  se_proj_relax <- sqrt(var(eic_proj_relax, na.rm = TRUE) / n)
  se_delta_reg <- sqrt(var(eic_delta_reg, na.rm = TRUE) / n)
  se_delta_relax <- sqrt(var(eic_delta_relax, na.rm = TRUE) / n)
  lower_proj_reg <- psi_reg - 1.96 * se_proj_reg
  lower_proj_relax <- psi_relax - 1.96 * se_proj_relax
  lower_delta_reg <- psi_reg - 1.96 * se_delta_reg
  lower_delta_relax <- psi_relax - 1.96 * se_delta_relax
  upper_proj_reg <- psi_reg + 1.96 * se_proj_reg
  upper_proj_relax <- psi_relax + 1.96 * se_proj_relax
  upper_delta_reg <- psi_reg + 1.96 * se_delta_reg
  upper_delta_relax <- psi_relax + 1.96 * se_delta_relax

  return(list(psi_reg = psi_reg,
              psi_relax = psi_relax,
              lower_proj_reg = lower_proj_reg,
              lower_proj_relax = lower_proj_relax,
              lower_delta_reg = lower_delta_reg,
              lower_delta_relax = lower_delta_relax,
              upper_proj_reg = upper_proj_reg,
              upper_proj_relax = upper_proj_relax,
              upper_delta_reg = upper_delta_reg,
              upper_delta_relax = upper_delta_relax))
}
