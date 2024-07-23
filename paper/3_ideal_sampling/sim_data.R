library(dplyr)
sim_data <- function(N_rwd, n_rct, n_rwd, g_rct, ate) {
  # EHR database ---------------------------------------------------------------
  UY <- rnorm(N_rwd, 0, 1)
  W1 <- rnorm(N_rwd, 0, 1)
  W2 <- rnorm(N_rwd, 0, 1)
  W3 <- rnorm(N_rwd, 0, 1)
  A <- rbinom(N_rwd, 1, plogis(-0.5+0.2*W1^2-0.6*W2^3))
  b <- 0.2+0.1*W1*(1-A)
  Y <- 2.5+0.9*W1+1.1*W2+2.7*W3+ate*A+UY+b
  ehr_data <- data.frame(W1, W2, W3, A, Y)

  # RWD 1: random sample from EHR database -------------------------------------
  rwd_rand_A1 <- slice_sample(ehr_data %>% filter(A == 1), n = (n_rwd/2), replace = FALSE)
  rwd_rand_A0 <- slice_sample(ehr_data %>% filter(A == 0), n = (n_rwd/2), replace = FALSE)
  rwd_rand_data <- rbind(rwd_rand_A1, rwd_rand_A0)
  rwd_rand_data$S <- 0

  # RWD 2: oracle sample -------------------------------------------------------
  rwd_oracle_A1 <- rwd_rand_A1
  rwd_oracle_A0 <- rwd_rand_A1; rwd_oracle_A0$A <- 0
  UY_oracle <- rnorm(n_rwd/2, 0, 1)
  b_oracle <- 0.2+0.1*rwd_oracle_A0$W1*(1-rwd_oracle_A0$A)
  rwd_oracle_A0$Y <- 2.5+0.9*rwd_oracle_A0$W1+1.1*rwd_oracle_A0$W2+2.7*rwd_oracle_A0$W3+ate*rwd_oracle_A0$A+UY_oracle+b_oracle
  rwd_oracle_data <- rbind(rwd_oracle_A1, rwd_oracle_A0)
  rwd_oracle_data$S <- 0

  # RCT ------------------------------------------------------------------------
  W1_rct <- rnorm(n_rct, 0, 1)
  W2_rct <- rnorm(n_rct, 0, 1)
  W3_rct <- rnorm(n_rct, 0, 1)
  A_rct <- rbinom(n_rct, 1, g_rct)
  UY_rct <- rnorm(n_rct, 0, 1)
  Y_rct <- 2.5+0.9*W1_rct+1.1*W2_rct+2.7*W3_rct+ate*A_rct+UY_rct
  rct_data <- data.frame(W1 = W1_rct,
                         W2 = W2_rct,
                         W3 = W3_rct,
                         A = A_rct,
                         Y = Y_rct)
  rct_data$S <- 1

  return(list(random = rbind(rwd_rand_data, rct_data),
              oracle = rbind(rwd_oracle_data, rct_data)))
}
