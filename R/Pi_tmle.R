#' TMLE update of the trial enrollment probabilities
#'
#' @export
#'
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param g A vector of estimated treatment probabilities,
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' @param tau A \code{list} containing the following elements:
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
#' @param Pi A \code{list} containing the following elements:
#' \item{pred}{The estimated trial enrollment probabilities;}
#' \item{A1}{The estimated trial enrollment probabilities under treatment;}
#' \item{A0}{The estimated trial enrollment probabilities under control.}
#' @param controls_only A logical indicating whether to learn only among
#' controls. This applies when the external data only has control-arm
#' observations.
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{Targeted estimates of trial enrollment probabilities;}
#' \item{A1}{Targeted estimates of trial enrollment probabilities under treatment;}
#' \item{A0}{Targeted estimates of trial enrollment probabilities under control.}
Pi_tmle <- function(S, W, A, g, tau, Pi, controls_only) {
  Pi_star <- Pi

  if (controls_only) {
    # clever covariate, controls only
    H0_n <- 1/(1-g[A == 0])*tau$A0[A == 0]

    # logistic submodel, controls only
    epsilon <- coef(glm(S[A == 0] ~ -1 + offset(qlogis(Pi$pred[A == 0])) + H0_n, family = "binomial"))
    epsilon[is.na(epsilon)] <- 0

    # TMLE update
    Pi_star$pred[A == 0] <- plogis(qlogis(Pi$pred[A == 0])+epsilon[1]*H0_n)
    Pi_star$A0[A == 0] <- plogis(qlogis(Pi$A0[A == 0])+epsilon[1]*H0_n)
  } else {
    # clever covariates, both treated and controls
    H1_n <- A/g*tau$A1
    H0_n <- (1-A)/(1-g)*tau$A0

    # logistic submodel, both treated and controls
    epsilon <- coef(glm(S ~ -1 + offset(qlogis(Pi$pred)) + H0_n + H1_n, family = "binomial"))
    epsilon[is.na(epsilon)] <- 0

    # TMLE updates
    Pi_star$pred <- plogis(qlogis(Pi$pred)+epsilon[1]*H0_n+epsilon[2]*H1_n)
    Pi_star$A0 <- plogis(qlogis(Pi$A0)+epsilon[1]*H0_n)
    Pi_star$A1 <- plogis(qlogis(Pi$A1)+epsilon[2]*H1_n)
  }

  return(Pi_star)
}
