sim_data <- function(n, controls_only, family, prop_miss) {
  UY <- rnorm(n, 0, 1)
  S <- rbinom(n, 1, 0.3)
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  A <- numeric(n)
  A[S == 1] <- rbinom(sum(S), 1, 0.5)

  if (controls_only) {
    A[S == 0] <- 0
  } else {
    A[S == 0] <- rbinom(n - sum(S), 1, plogis(1.2*W1[S==0]-0.9*W2[S==0]))
  }

  if (family == "gaussian") {
    Y <- -0.5-0.8*W1-1.1*W2+1.5*A+UY+(1-S)*(0.9+2.6*W1)
  } else if (family == "binomial") {
    Y <- rbinom(n, 1, plogis(-0.9*W1-0.5*W2+0.5*A+(1-S)*(0.9-0.3*W1)))
  }

  if (prop_miss > 0) {
    Y[runif(n) < prop_miss] <- NA
  }
  return(data.frame(S = S,
                    W1 = W1,
                    W2 = W2,
                    A = A,
                    Y = Y))
}

# This function is adapted from CRAN package `rlist`
list.expand <- function(...) {
  args <- list(...)
  expand_args <- lapply(args, seq_along)
  expand_df <- do.call(expand.grid, expand_args)
  .mapply(function(...) {
    mapply(`[[`, args, list(...), USE.NAMES = TRUE, SIMPLIFY = FALSE)
  }, expand_df, NULL)
}
