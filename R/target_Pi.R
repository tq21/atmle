#' @title TMLE update of the trial enrollment probabilities
#'
#' @description This function updates the trial enrollment probabilities
#' using the targeted maximum likelihood estimation (TMLE).
#'
#' @keywords TMLE
#'
#' @param S A vector of study indicators, \eqn{S=1} for RCT, \eqn{S=0} for RWD.
#' @param W A matrix of baseline covariates.
#' @param A A vector of treatment indicators, \eqn{A=1} for treatment-arm,
#' \eqn{A=0} for control-arm.
#' @param g A vector of estimated treatment probabilities,
#' \eqn{g(A\mid W)=\mathbb{P}(A=1\mid W)}.
#' @param tau Returns from the \code{\link{learn_tau}} function.
#' @param Pi Returns from the \code{\link{learn_Pi}} function.
#' @param controls_only A logical indicating whether the external data has only
#' control-arm observations.
#' @param target_gwt If \code{TRUE}, the treatment mechanism is moved from the
#' denominator of the clever covariate to the weight when fitting the TMLE
#' submodel.
#' @param Pi_bounds A numeric vector of lower and upper bounds for the
#' trial enrollment probabilities. The first element is the lower bound,
#' and the second element is the upper bound.
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{Targeted estimates of trial enrollment probabilities;}
#' \item{A1}{Targeted estimates of trial enrollment probabilities under treatment;}
#' \item{A0}{Targeted estimates of trial enrollment probabilities under control.}
target_Pi <- function(S,
                      W,
                      A,
                      Y,
                      g1W,
                      Pi,
                      theta_WA,
                      tau_S,
                      controls_only,
                      target_gwt,
                      weights,
                      Pi_bounds) {

  Pi_star <- Pi

  if (controls_only) {
    if (target_gwt) {
      wt <- 1/(1-g1W[A == 0])
      H0_n <- tau_S$cate_W0[A == 0]
    } else {
      wt <- rep(1, length(A[A == 0]))
      H0_n <- 1/(1-g1W[A == 0])*tau_S$cate_W0[A == 0]
    }

    # logistic submodel, controls only
    epsilon <- coef(glm(S[A == 0] ~ -1 + offset(qlogis(Pi$pred[A == 0])) + H0_n,
                        family = "quasibinomial", weights = wt
    ))
    epsilon[is.na(epsilon)] <- 0

    # TMLE update
    if (target_gwt) {
      Pi_star$pred[A == 0] <- .bound(plogis(qlogis(Pi$pred[A == 0]) + epsilon[1] * H0_n), Pi_bounds)
      Pi_star$A0[A == 0] <- .bound(plogis(qlogis(Pi$A0[A == 0]) + epsilon[1] * H0_n), Pi_bounds)
    } else {
      Pi_star$pred[A == 0] <- .bound(plogis(qlogis(Pi$pred[A == 0]) + epsilon[1] * H0_n), Pi_bounds)
      Pi_star$A0[A == 0] <- .bound(plogis(qlogis(Pi$A0[A == 0]) + epsilon[1] * H0_n), Pi_bounds)
    }
  } else {
    if (target_gwt) {
      wt <- A/g1W+(1-A)/(1-g1W)
      H1_n <- tau_S$cate_W1*A
      H0_n <- tau_S$cate_W0*(1-A)
    } else {
      wt <- rep(1, length(A))
      H1_n <- A/g1W*tau_S$cate_W1
      H0_n <- (1-A)/(1-g1W)*tau_S$cate_W0
    }

    # logistic submodel, both treated and controls
    epsilon <- coef(glm(S ~ -1 + offset(qlogis(Pi$pred)) + H0_n + H1_n,
                        family = "quasibinomial", weights = wt
    ))
    epsilon[is.na(epsilon)] <- 0

    # TMLE updates
    if (target_gwt) {
      Pi_star$pred <- .bound(plogis(qlogis(Pi$pred) + epsilon[1] * H0_n + epsilon[2] * H1_n), Pi_bounds)
      Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * tau_S$cate_W0), Pi_bounds)
      Pi_star$A1 <- .bound(plogis(qlogis(Pi$A1) + epsilon[2] * tau_S$cate_W1), Pi_bounds)
    } else {
      Pi_star$pred <- .bound(plogis(qlogis(Pi$pred) + epsilon[1] * H0_n + epsilon[2] * H1_n), Pi_bounds)
      Pi_star$A0 <- .bound(plogis(qlogis(Pi$A0) + epsilon[1] * tau_S$cate_W0/(1-g1W)), Pi_bounds)
      Pi_star$A1 <- .bound(plogis(qlogis(Pi$A1) + epsilon[2] * tau_S$cate_W1/g1W), Pi_bounds)
    }
  }

  # update relevant parts of tau_S
  tau_S$pseudo_outcome <- ifelse(abs(S-Pi_star$pred) < 1e-10, 0, (Y-theta_WA)/(S-Pi_star$pred))
  tau_S$pseudo_weights <- (S-Pi_star$pred)^2*weights

  return(list(Pi = Pi_star,
              tau_S = tau_S))
}
