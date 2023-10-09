#' @description Function to compute bootstrap confidence interval for
#' the bias parameter
bootstrap_psi_pound <- function(tau, W, Pi, B = 1000) {
  psi_pound_var <- var(unlist(map(seq(B), function(b) {
    # sample rows
    psi_pound_samp_idx <- sample(nrow(tau$x_basis), replace = TRUE)

    # predictions
    A1 <- as.numeric(tau$x_basis_A1[psi_pound_samp_idx,] %*% matrix(tau$coefs))
    A0 <- as.numeric(tau$x_basis_A0[psi_pound_samp_idx,] %*% matrix(tau$coefs))

    return(mean((1-Pi$A0[psi_pound_samp_idx])*A0-(1-Pi$A1[psi_pound_samp_idx])*A1))
  })))

  return(sqrt(psi_pound_var))
}

#' @description Function to compute bootstrap confidence interval for
#' the pooed-ATE parameter, does not apply when tmle is used to estimate the
#' pooled-ATE parameter
bootstrap_psi_tilde <- function(W, psi_tilde, B = 1000) {
  psi_tilde_var <- var(unlist(map(seq(B), function(b) {
    # sample rows
    psi_tilde_samp_idx <- sample(nrow(W), replace = TRUE)

    # predict
    return(mean(as.numeric(as.matrix(cbind(1, W)) %*% matrix(psi_tilde$coefs))))
  })))

  return(sqrt(psi_tilde_var))
}
