#' @title Function to compute bootstrap confidence interval for the bias
#' estimate.
#'
#' @keywords internal
#'
#' @importFrom purrr map
#'
#' @param A \code{list} containing the following elements:
#' \item{A1}{A numeric vector of the estimated counterfactual conditional
#' effects under treatment;}
#' \item{A0}{A numeric vector of the estimated counterfactual conditional
#' effects under control;}
#' \item{x_basis}{A numeric matrix of the working model bases;}
#' \item{x_basis_A1}{A numeric matrix of the counterfactual working model bases
#' under treatment;}
#' \item{x_basis_A0}{A numeric matrix of the counterfactual working model bases
#' under control;}
#' \item{pred}{A numeric vector of estimated conditional effects;}
#' \item{coefs}{A numeric vector of the working model coefficients.}
#' @param W A matrix of baseline covariates.
#' @param Pi A \code{list} containing the following elements:
#' \item{pred}{The estimated trial enrollment probabilities;}
#' \item{A1}{The estimated trial enrollment probabilities under treatment;}
#' \item{A0}{The estimated trial enrollment probabilities under control.}
#' @param B The number of bootstrap samples.
bootstrap_psi_pound <- function(tau, W, Pi, B = 1000) {
  psi_pound_var <- var(unlist(map(seq(B), function(b) {
    # sample rows
    psi_pound_samp_idx <- sample(nrow(tau$x_basis), replace = TRUE)

    # predictions
    A1 <- as.numeric(tau$x_basis_A1[psi_pound_samp_idx, ] %*% matrix(tau$coefs))
    A0 <- as.numeric(tau$x_basis_A0[psi_pound_samp_idx, ] %*% matrix(tau$coefs))

    return(mean((1 - Pi$A0[psi_pound_samp_idx]) * A0 - (1 - Pi$A1[psi_pound_samp_idx]) * A1))
  })))

  return(sqrt(psi_pound_var))
}

#' @title Function to compute bootstrap confidence interval for
#' the pooled-ATE estimate.
#'
#' @importFrom purrr map
#'
#' @param W A matrix of baseline covariates.
#' @param psi_tilde A \code{list} containing the following elements:
#' \item{pred}{A numeric vector of the estimated conditional effects;}
#' \item{x_basis}{A numeric matrix of the working model bases;}
#' \item{coefs}{A numeric vector of the working model coefficients.}
#' @param B The number of bootstrap samples.
bootstrap_psi_tilde <- function(W, psi_tilde, B = 1000) {
  psi_tilde_var <- var(unlist(map(seq(B), function(b) {
    # sample rows
    psi_tilde_samp_idx <- sample(nrow(W), replace = TRUE)

    # predict
    return(mean(as.numeric(as.matrix(cbind(1, W)) %*% matrix(psi_tilde$coefs))))
  })))

  return(sqrt(psi_tilde_var))
}
