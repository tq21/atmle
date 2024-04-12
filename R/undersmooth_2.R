undersmooth_2 <- function(S, W, A, X, Y,
                          Pi, theta, g,
                          X_A1, X_A0,
                          pseudo_weights,
                          X_unpenalized,
                          controls_only,
                          v_folds,
                          target_gwt,
                          Pi_bounds) {

  # fit regular HAL
  fit <- fit_relaxed_hal(X = X, Y = Y,
                         X_unpenalized = X_unpenalized,
                         X_weak_penalized = NULL,
                         X_weak_penalized_level = 0,
                         family = "gaussian",
                         weights = pseudo_weights,
                         relaxed = FALSE,
                         v_folds = v_folds,
                         hal_args = list())
  hal_fit <- fit$hal_fit

  # undersmooth to get larger model, but still solve PnD*
  n <- length(S)
  lasso_obj <- hal_fit$lasso_fit
  lambda_min <- lasso_obj$lambda.min
  lambda_seq <- lasso_obj$lambda[lasso_obj$lambda <= lambda_min]

  undersmooth_obj <- map(lambda_seq, function(lambda) {
    # get selected basis under the current lambda
    coefs <- as.numeric(coef(lasso_obj, s = lambda))
    selected_idx <- which(coefs != 0)
    x_basis <- cbind(1, hal_fit$x_basis)[, selected_idx, drop = FALSE]

    # get penalized fit
    pred <- predict(hal_fit,
                    new_data = X,
                    new_X_unpenalized = X_unpenalized,
                    type = "response",
                    s = lambda)

    # design matrices
    hal_basis_list <- hal_fit$basis_list[coefs[-1] != 0]
    if (controls_only) {
      x_basis_A1 <- x_basis_A0 <- make_counter_design_matrix(hal_basis_list,
                                                             as.matrix(W),
                                                             X_unpenalized = X_unpenalized)
    } else {
      x_basis_A1 <- make_counter_design_matrix(hal_basis_list,
                                               as.matrix(cbind(W, A = 1)),
                                               X_unpenalized = X_unpenalized)
      x_basis_A0 <- make_counter_design_matrix(hal_basis_list,
                                               as.matrix(cbind(W, A = 0)),
                                               X_unpenalized = X_unpenalized)
    }

    # predictions
    A1 <- predict(hal_fit,
                  new_data = X_A1,
                  new_X_unpenalized = X_unpenalized,
                  type = "response",
                  s = lambda)
    A0 <- predict(hal_fit,
                  new_data = X_A0,
                  new_X_unpenalized = X_unpenalized,
                  type = "response",
                  s = lambda)

    # make tau object
    tau <- list(A1 = A1,
                A0 = A0,
                x_basis = as.matrix(x_basis),
                x_basis_A1 = x_basis_A1,
                x_basis_A0 = x_basis_A0,
                pred = pred,
                coefs = beta,
                non_zero = selected_idx,
                pseudo_outcome = Y,
                pseudo_weights = pseudo_weights)

    # compute psi_pound_est
    if (controls_only) {
      psi_pound_est <- mean((1-Pi$A0)*tau$A0)
    } else {
      psi_pound_est <- mean((1-Pi$A0)*tau$A0-(1-Pi$A1)*tau$A1)
    }

    # target Pi
    Pi <- Pi_tmle(S = S, W = W, A = A,
                  g = g, tau = tau, Pi = Pi,
                  controls_only = controls_only,
                  target_gwt = target_gwt,
                  Pi_bounds = Pi_bounds)

    # evaluate EIC
    tryCatch({
      D_star <- get_eic_psi_pound(Pi = Pi,
                                  tau = tau,
                                  g = g,
                                  theta = theta,
                                  psi_pound_est = psi_pound_est,
                                  S = S,
                                  A = A,
                                  Y = Y,
                                  n = n,
                                  controls_only = controls_only)
      PnD_star <- mean(D_star)
      tol <- sqrt(mean(D_star^2))/(log(n)*sqrt(n))
      return(list(PnD_star = PnD_star,
                  tol = tol,
                  tau = tau))
    }, error = function(e) {
      return(list(PnD_star = Inf,
                  tol = Inf,
                  tau = tau))
    })
  })

  # find the first lambda that gives PnD* <= tol
  PnD_star_all <- map_vec(undersmooth_obj, function(.x) .x$PnD_star)
  tol_all <- map_vec(undersmooth_obj, function(.x) .x$tol)
  idxes <- which(abs(PnD_star_all) <= tol_all)
  if (length(idxes) == 0) {
    idx <- which.min(abs(PnD_star_all)) # 1
  } else {
    idx <- idxes[1]
  }

  return(undersmooth_obj[[idx]]$tau)
}
