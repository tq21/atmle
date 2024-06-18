sim_data <- function(n, A_counter = NULL) {

  W1 <- round(runif(n, 2, 6), 2)
  W2 <- round(rnorm(n, 10, sqrt(10)), 2)

  if (!is.null(A_counter)) {
    A <- rep(A_counter, n)
  } else {
    A <- rbinom(n, 1, 0.5)
  }

  # failure hazard
  lambda_fun <- function(t, A, W1, W2) {
    if (t == 9) {
      return(1)
    } else if (t %in% 1:8) {
      return(plogis(-8-0.75*A+0.3*W1+0.25*W2))
    }
  }
  lambda_fun <- Vectorize(lambda_fun)

  # censoring hazard
  lambda_c_fun <- function(t, A, W1, W2) {
    if (t == 1) {
      return(0)
    } else if (W1 > 4.5 & A == 1) {
      return(0.25)
    } else if (W1 > 3.5 & W1 <= 4.5 & A == 1) {
      return(0.2)
    } else if (W1 > 2.5 & W1 <= 3.5 & A == 1) {
      return(0.05)
    } else if (W1 > 3.5 & A == 0) {
      return(0)
    } else if (W1 > 2.5 & W1 <= 3.5 & A == 0) {
      return(0.25)
    } else {
      return(0.05)
    }
  }
  #lambda_c_fun <- function(t, A, W1, W2) return(0)
  lambda_c_fun <- Vectorize(lambda_c_fun)

  dt <- data.table(id = rep(seq(n), each = 9),
                   t = rep(1:9, n),
                   W1 = rep(W1, each = 9),
                   W2 = rep(W2, each = 9),
                   A = rep(A, each = 9))
  dt[, `:=` (lambda = lambda_fun(t, A, W1, W2),
             lambda_c = lambda_c_fun(t, A, W1, W2))]
  dt[, `:=` (surv = cumprod(1 - lambda),
             surv_c = cumprod(1 - lambda_c)), by = id]
  dt[, `:=` (prob = lambda * shift(surv, fill = 1),
             prob_c = lambda_c * shift(surv_c, fill = 1)), by = id]
  dt[, `:=` (T = rbinom(.N, 1, prob),
             C = rbinom(.N, 1, prob_c))]
  dt[, `:=` (T_t = T * t,
             C_t = C * t)]
  dt[T_t == 0, T_t := 999]
  dt[C_t == 0, C_t := 9]
  dt[, T_t := min(T_t), by = id]
  dt[, C_t := min(C_t), by = id]
  dt[, T_tilde := min(T_t, C_t), by = id]
  dt[, Delta := as.numeric(T_t <= C_t), by = id]

  return(as.data.frame(dt[t == 1, .(W1, W2, A, T_tilde, Delta)]))
}

library(devtools)
library(data.table)
load_all()
library(hal9001)
library(glmnet)
library(purrr)

t0 <- 6
data_A1 <- sim_data(10000, A_counter = 1)
data_A0 <- sim_data(10000, A_counter = 0)
mean(data_A1$T_tilde >= t0)-mean(data_A0$T_tilde >= t0)

n <- 2000
data <- sim_data(n)
res <- atmle_surv(data = data,
                  S = "T_tilde",
                  W = c("W1", "W2"),
                  A = "A",
                  T_tilde = "T_tilde",
                  Delta = "Delta",
                  t0 = t0,
                  g_rct = 0.5,
                  controls_only = FALSE,
                  g_method = "glm",
                  G_bar_method = "glm",
                  lambda_method = "glm",
                  theta_method = "glm")
