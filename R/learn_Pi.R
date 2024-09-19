#' @title Learn nuisance function: trial enrollment mechanism
#'
#' @description Function to learn the trial enrollment mechanism
#' \eqn{\Pi(S\mid W,A)=\mathbb{P}(S=1\mid W,A)}.
#'
#' @keywords nuisance
#'
#' @param g Returned list from `learn_g()` function.
#'
#' @returns A \code{list} containing the following elements:
#' \item{pred}{The estimated trial enrollment probabilities;}
#' \item{A1}{The estimated trial enrollment probabilities under treatment;}
#' \item{A0}{The estimated trial enrollment probabilities under control.}
learn_Pi <- function(g,
                     A,
                     Pi_bounds) {
  # P(S=1|W,A=1)
  Pi_A1 <- (g$pred_A_S1_W * g$pred_S_W) / g$pred

  # P(S=1|W,A=0)
  Pi_A0 <- ((1 - g$pred_A_S1_W) * g$pred_S_W) / (1 - g$pred)

  # P(S|W,A)
  Pi <- numeric(length(A)); Pi[A == 1] <- Pi_A1[A == 1]; Pi[A == 0] <- Pi_A0[A == 0]

  return(list(pred = .bound(Pi, Pi_bounds),
              A1 = .bound(Pi_A1, Pi_bounds),
              A0 = .bound(Pi_A0, Pi_bounds)))
}
