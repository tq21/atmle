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
bootstrap_ci <- function(tau,
                         T_working,
                         Pi,
                         tau_weights,
                         B = 500) {

  boot_ests <- map_vec(seq(B), function(b) {
    # sample rows
    samp_idx <- sample(nrow(tau$x_basis), replace = TRUE)

    # learn psi pound
    fit_tau <- glm.fit(x = tau$x_basis[samp_idx, ],
                       y = tau$pseudo_outcome[samp_idx],
                       family = gaussian(),
                       weights = tau$pseudo_weights,
                       intercept = FALSE)
    beta_tau <- coef(fit_tau); beta_tau[is.na(beta_tau)] <- 0
    A1 <- as.numeric(tau$x_basis_A1[samp_idx, ] %*% matrix(beta_tau))
    A0 <- as.numeric(tau$x_basis_A0[samp_idx, ] %*% matrix(beta_tau))

    # learn psi tilde
    fit_T_working <- glm.fit(x = T_working$x_basis[samp_idx, ],
                             y = T_working$pseudo_outcome[samp_idx],
                             family = gaussian(),
                             weights = T_working$pseudo_weights,
                             intercept = FALSE)
    beta_T_working <- coef(fit_T_working); beta_T_working[is.na(beta_T_working)] <- 0
    pred <- as.numeric(T_working$x_basis[samp_idx, ] %*% matrix(beta_T_working))

    # get point estimate
    if (controls_only) {
      psi_pound_est <- mean((1-A0)*A0)
    } else {
      psi_pound_est <- mean((1-Pi$A0[samp_idx])*A0-(1-Pi$A1[samp_idx])*A1)
    }
    psi_tilde_est <- mean(pred)
    psi <- psi_tilde_est - psi_pound_est

    return(psi)
  })

  return(list(lower = as.numeric(quantile(boot_ests, probs = 0.025)),
              upper = as.numeric(quantile(boot_ests, probs = 0.975))))
}
