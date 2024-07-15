library(data.table)
sim_data <- function(n, A_counter = NULL) {

  W1 <- round(runif(n, 0, 4), 2)
  W2 <- round(runif(n, 0, 4), 2)

  if (!is.null(A_counter)) {
    A <- rep(A_counter, n)
  } else {
    A <- rbinom(n, 1, 0.5)
  }

  # failure hazard
  lambda_fun <- function(t, A, W1, W2) {
    return(plogis(-0.5*A-0.5*W1))
  }
  lambda_fun <- Vectorize(lambda_fun)

  # censoring hazard
  lambda_c_fun <- function(t, A, W1, W2) {
    if (t == 1) {
      return(0)
    } else {
      return(0.1)
    }
  }
  if (!is.null(A_counter)) {
    lambda_c_fun <- function(t, A, W1, W2) return(0)
  }
  lambda_c_fun <- Vectorize(lambda_c_fun)

  dt <- data.table(id = rep(seq(n), each = 5),
                   t = rep(1:5, n),
                   W1 = rep(W1, each = 5),
                   W2 = rep(W2, each = 5),
                   A = rep(A, each = 5))
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
  dt[T_t == 0, T_t := 6]
  dt[C_t == 0, C_t := 5]
  dt[, T_t := min(T_t), by = id]
  dt[, C_t := min(C_t), by = id]
  dt[, T_tilde := min(T_t, C_t), by = id]
  dt[, Delta := as.numeric(T_t <= C_t), by = id]

  return(as.data.frame(dt[t == 1, .(W1, W2, A, T_tilde, Delta)]))
}

get_true_cate <- function(data, t0) {
  W1 <- data$W1; W2 <- data$W2
  cate <- (1-plogis(-0.5-0.5*W1))^t0 - (1-plogis(-0.5*W1))^t0
  return(cate)
}

data <- sim_data(2000)
summary(as.factor(data$T_tilde))
summary(as.factor(data$Delta))
