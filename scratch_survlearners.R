library(survlearners)
set.seed(123)
t0 <- 3

n <- 1000
B <- 200

res <- numeric(B)

for (b in 1:B) {
  print(b)
  data <- sim_data(n)
  fit <- surv_rl_lasso(X = as.matrix(data[, c("W1", "W2")]),
                       Y = data$T_tilde,
                       W = data$A,
                       D = data$Delta,
                       t0 = t0,
                       W.hat = rep(0.5, n),
                       cen.fit = "survival.forest")
  res[b] <- mean(mean(as.numeric(predict(fit))))
}

hist(res, main = "IPCW R-learner")
abline(v = truth, col = "red")
