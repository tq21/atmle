library(torch)
atmle_torch <- function(data,
                        S,
                        W,
                        A,
                        Y,
                        controls_only,
                        family,
                        lr = 1e-3,
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
                        cross_fit_nuisance = TRUE,
                        verbose = TRUE,
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

  # learn nuisance parts
  if (verbose) cat("learning \U03B8(W,A)=E(Y|W,A)...")
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
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning \U03B8\U0303(W)=E(Y|W)...")
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
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning g(A=1|W)=P(A=1|W)...")
  g <- learn_g(
    S = S,
    W = W,
    A = A,
    method = g_method,
    controls_only = controls_only,
    v_folds = v_folds,
    g_bounds = g_bounds
  )
  if (verbose) cat("Done!\n")

  if (verbose) cat("learning \U03A0(S=1|W,A)=P(S=1|W,A)...")
  Pi <- learn_Pi(
    g = g,
    A = A,
    Pi_bounds = Pi_bounds
  )
  if (verbose) cat("Done!\n")

  if (sum(delta) < n) {
    # outcome has missing
    if (verbose) cat("learning g(\U0394=1|S,W,A)=P(\U0394=1|S,W,A)...")
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

    if (verbose) cat("learning g(\U0394=1|W,A)=P(\U0394=1|W,A)...")
    g_delta_tilde <- learn_g_delta_tilde(
      W = W,
      A = A,
      delta = delta,
      method = g_delta_method,
      folds = folds,
      g_bounds = g_bounds
    )
    if (verbose) cat("Done!\n")
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

  browser()
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
    patience = patience
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
    patience = patience
  )
}
