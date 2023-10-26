sim_four_covs <- function(ate, n_rct, n_rwd, g_rct, bias, controls_only) {
  n <- n_rct + n_rwd

  # simulate trial data --------------------------------------------------------
  # error
  UY <- rnorm(n_rct, 0, 1.5)

  # baseline covariates
  W1 <- rnorm(n_rct, 0, 1)
  W2 <- rnorm(n_rct, 0, 1)
  W3 <- rnorm(n_rct, 0, 1)
  W4 <- rnorm(n_rct, 0, 1)

  # treatment
  A <- rbinom(n_rct, 1, g_rct)

  # outcome
  Y <- -5.2+3.1*W1+1.1*W2+0.9*W3+1.3*W4+ate*A+UY

  # trial data frame
  rct_dt <- data.frame(S = 1,
                       W1 = W1,
                       W2 = W2,
                       W3 = W3,
                       W4 = W4,
                       A = A,
                       Y = Y)

  # simulate external data -----------------------------------------------------
  large_n <- 10*n_rwd

  # error
  U_bias <- rnorm(large_n, 0, 0.01)
  UY <- rnorm(large_n, 0, 1.5)

  # baseline covariates
  W1 <- rnorm(large_n, 0, 1)
  W2 <- rnorm(large_n, 0, 1)
  W3 <- rnorm(large_n, 0, 1)
  W4 <- rnorm(large_n, 0, 1)

  # treatment
  A <- rbinom(large_n, 1, plogis(0.9*W1+0.4*W2))

  # bias term
  b <- NULL
  if (is.numeric(bias)) {
    b <- bias
  } else if (bias == "param_simple") {
    b <- 2.6*W1
  } else if (bias == "param_complex") {
    b <- 2.4+5.3*W1+2.7*W2
  } else if (bias == "HAL") {
    b <- 0.3+8*W1^2.3*W2
  } else if (bias == "sinusoidal") {
    b <- 3.8*W2*as.numeric(W1 < 0.5)*sin(pi/2*abs(W1))+
      4*as.numeric(W2 > 0.7)*cos(pi/2*abs(W2))
  }

  # outcome
  Y <- -5.2+3.1*W1+1.1*W2+0.9*W3+1.3*W4+ate*A+UY+b+U_bias

  # external data frame
  rwd_dt <- data.frame(S = 0,
                       W1 = W1,
                       W2 = W2,
                       W3 = W3,
                       W4 = W4,
                       A = A,
                       Y = Y)

  # sample external data
  rwd_dt_samp <- NULL
  if (controls_only) {
    # only sample external controls
    rwd_dt_samp <- rwd_dt[sample(which(rwd_dt$A == 0), n_rwd), ]
  } else {
    # sample external data, half treated and half control
    #rwd_dt_trt <- rwd_dt[sample(which(rwd_dt$A == 1), n_rwd/2), ]
    #rwd_dt_con <- rwd_dt[sample(which(rwd_dt$A == 0), n_rwd/2), ]
    #rwd_dt_samp <- rbind(rwd_dt_trt, rwd_dt_con)
    rwd_dt_samp <-rwd_dt[sample(nrow(rwd_dt), n_rwd), ]
  }

  # data frames combining RCT and RWD
  data <- rbind(rct_dt, rwd_dt_samp)

  return(data)
}
