naive_atmle <- function(data,
                        S,
                        W,
                        A,
                        Y,
                        family,
                        g_method = "glmnet",
                        theta_tilde_method = "glmnet",
                        cross_fit_nuisance = TRUE,
                        v_folds = 5,
                        g_bounds = c(0.01, 0.99),
                        theta_bounds = c(-Inf, Inf)) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # extract variables ----------------------------------------------------------
  S <- data[[S]]
  W <- data[, W, drop = FALSE]
  A <- data[[A]]
  Y <- data[[Y]]
  n <- sum(S)

  # cv.glmnet only works when design matrix has at least 2 columns
  # append dummy column of ones if necessary
  if (ncol(W) == 1) {
    W <- cbind(W, 1)
  }

  # cross fitting schemes
  if (family == "gaussian") {
    cv_strata <- A[S == 1]
  } else if (family == "binomial") {
    cv_strata <- paste0(A[S == 1], "-", Y[S == 1])
  }

  suppressWarnings({
    folds <- make_folds(
      n = n, V = v_folds,
      strata_ids = as.integer(factor(cv_strata))
    )
  })

  # estimate nuisance: g
  g <- learn_g(
    W = W[S == 1, ],
    A = A[S == 1],
    method = g_method,
    folds = folds,
    g_bounds = g_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  # estimate nuisance: theta_tilde
  theta_tilde <- learn_theta_tilde(
    W = W[S == 1, ],
    Y = Y[S == 1],
    delta = rep(1, n),
    method = theta_tilde_method,
    folds = folds,
    family = family,
    theta_bounds = theta_bounds,
    cross_fit_nuisance = cross_fit_nuisance
  )

  # learn working model
  # R-transformations
  pseudo_outcome <- ifelse(abs(A[S == 1] - g) < 1e-10, 0, (Y[S == 1] - theta_tilde) / (A[S == 1] - g))
  pseudo_weights <- (A[S == 1] - g)^2

  fit <- cv.glmnet(x = as.matrix(W[S == 1, ]), y = pseudo_outcome, weights = pseudo_weights,
                   family = "gaussian", keep = TRUE, nfolds = v_folds, alpha = 1, relax = TRUE)
  non_zero <- which(as.numeric(coef(fit, s = "lambda.min", gamma = 0)) != 0)
  coefs <- coef(fit, s = "lambda.min", gamma = 0)[non_zero]
  x_basis <- as.matrix(cbind(1, W)[, non_zero, drop = FALSE])
  pred <- as.numeric(x_basis %*% matrix(coefs))

  return(mean(pred))
}
