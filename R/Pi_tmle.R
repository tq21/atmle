# function to perform TMLE update of Pi
Pi_tmle <- function(S, W, A, g, tau, Pi) {
  H1_n <- A/g*tau$A1
  H0_n <- (1-A)/(1-g)*tau$A0

  # logistic submodel
  suppressWarnings(
    epsilon <- coef(glm(S ~ -1 + offset(qlogis(Pi$pred)) + H0_n + H1_n, family = "binomial"))
  )
  epsilon[is.na(epsilon)] <- 0

  # TMLE updates
  Pi_star <- list()
  Pi_star$pred <- plogis(qlogis(Pi$pred)+epsilon[1]*H0_n+epsilon[2]*H1_n)
  Pi_star$A0 <- plogis(qlogis(Pi$A0)+epsilon[1]*H0_n)
  Pi_star$A1 <- plogis(qlogis(Pi$A1)+epsilon[2]*H1_n)

  return(Pi_star)
}
