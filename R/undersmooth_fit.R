undersmooth_fit <- function(hal_fit,
                            S, W, A, X, Y,
                            Pi, theta,
                            X_unpenalized,
                            controls_only,
                            relaxed) {
  n <- length(S)
  lasso_obj <- hal_fit$lasso_fit
  lambda_min <- lasso_obj$lambda.min
  lambda_seq <- lasso_obj$lambda[lasso_obj$lambda <= lambda_min]

  i <- 1
  PnD_star <- 1
  tol <- 0
  while (PnD_star > tol && i <= length(lambda_seq)) {
    if (i > 1) {
      tau_good <- tau
    }

    # current lambda
    lambda <- lambda_seq[i]

    # get selected basis under the current lambda
    coefs <- as.numeric(coef(lasso_obj, s = lambda))
    selected_idx <- which(coefs != 0)
    x_basis <- cbind(1, hal_fit$x_basis)[, selected_idx, drop = FALSE]

    if (relaxed) {
      # get relaxed fit
      fit <- glm.fit(x = x_basis, y = Y, family = gaussian(), intercept = FALSE)
      beta <- coef(fit)
      beta[is.na(beta)] <- 0
      pred <- as.vector(x_basis %*% beta)
    } else {
      # get penalized fit
      coefs <- as.numeric(coef(lasso_obj, s = lambda))
      beta <- coefs[coefs != 0]
      pred <- predict(hal_fit,
        new_data = X,
        new_X_unpenalized = X_unpenalized,
        type = "response",
        s = lambda
      )
    }

    # design matrices
    hal_basis_list <- hal_fit$basis_list[coefs[-1] != 0]
    x_basis_A1 <- make_counter_design_matrix(hal_basis_list,
      as.matrix(cbind(W, A = 1)),
      X_unpenalized = X_unpenalized
    )
    x_basis_A0 <- make_counter_design_matrix(hal_basis_list,
      as.matrix(cbind(W, A = 0)),
      X_unpenalized = X_unpenalized
    )

    # predictions
    A1 <- as.numeric(x_basis_A1 %*% matrix(beta))
    A0 <- as.numeric(x_basis_A0 %*% matrix(beta))

    # make tau object
    tau <- list(
      A1 = A1,
      A0 = A0,
      x_basis = as.matrix(x_basis),
      x_basis_A1 = x_basis_A1,
      x_basis_A0 = x_basis_A0,
      pred = pred,
      coefs = beta,
      non_zero = selected_idx
    )

    if (i == 1) {
      tau_good <- tau
    }

    # compute psi_pound_est
    if (controls_only) {
      psi_pound_est <- mean((1 - Pi$A0) * tau$A0)
    } else {
      psi_pound_est <- mean((1 - Pi$A0) * tau$A0 - (1 - Pi$A1) * tau$A1)
    }

    # evaluate EIC
    IM <- t(tau$x_basis) %*% diag((Pi$pred * (1 - Pi$pred))) %*% tau$x_basis / n
    D <- tau$x_basis %*% solve(IM) * (S - Pi$pred) * (Y - theta - (S - Pi$pred) * tau$pred)
    D_star <- NULL
    if (ncol(D) > 1) {
      D_star <- rowSums(D %*% diag(colMeans((1 - Pi$A0) * tau$x_basis_A0))) - rowSums(D %*% diag(colMeans((1 - Pi$A1) * tau$x_basis_A1)))
    } else {
      D_star <- rowSums(D * colMeans((1 - Pi$A0) * tau$x_basis_A0)) - rowSums(D * colMeans((1 - Pi$A1) * tau$x_basis_A1))
    }

    PnD_star <- mean(D_star)
    tol <- sqrt(mean(D_star^2)) / (log(n) * sqrt(n))
    i <- i + 1
  }

  print("undersmoothed: " %+% i)
  print("number of non-zero coefficients: " %+% length(tau$non_zero))

  return(tau)
}
