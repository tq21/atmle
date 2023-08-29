library(ggplot2)

load_all()

generate_data <- function(N, bA, S=NULL, A=NULL){
  # U
  UY <- rnorm(N, 0, 1)

  # W
  W1 <- rnorm(N)
  W2 <- rnorm(N)
  W3 <- rnorm(N)
  W4 <- rnorm(N)

  # A
  if (is.null(A)) {
    # A <- rbinom(N, 1, plogis(-0.1+0.2*W1+0.5*W2-0.1*W3))
    A <- rbinom(N, 1, 0.5)
  }

  # S
  if (is.null(S)) {
    S <- rbinom(N, 1, 0.5)
  }

  # Y
  Y <- 0.3+10*S*A+bA*A+0.5*W1+0.3*W3-0.5*W4+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

n <- 1000
bA <- 1.3
S1_true_ate_A <- mean(generate_data(100000, bA, rep(1, 100000), rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, rep(1, 100000), rep(0, 100000))$Y)
S0_true_ate_A <- mean(generate_data(100000, bA, rep(0, 100000), rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, rep(0, 100000), rep(0, 100000))$Y)
truth <- 0.5*S0_true_ate_A-0.5*S1_true_ate_A
truth <- -5

test <- function(n, bA) {
  data <- generate_data(n, bA)
  S_node <- 1
  W_node <- c(2, 3, 4, 5)
  A_node <- 6
  Y_node <- 7
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]

  # learn nuisance parts
  theta <- learn_theta(W, A, Y, method = "lasso")
  Pi <- learn_Pi(S, W, A, method = "lasso")
  g <- learn_g_tmp(W, A, method = "lasso")

  # learn initial estimate of working model tau
  tau <- learn_tau(S, W, A, Y, Pi$pred, theta, method = "lasso")

  for (i in 1:1) {
    # TMLE to target Pi
    Pi_star <- Pi_tmle(S, W, A, g, tau, Pi)

    # re-learn working model tau with targeted Pi
    tau_star <- learn_tau(S, W, A, Y, Pi_star$pred, theta, method = "lasso")

    Pi <- Pi_star
    tau <- tau_star
  }

  psi_pound_est <- mean((1-Pi_star$A0)*tau_star$A0-(1-Pi_star$A1)*tau_star$A1)
  psi_pound_eic <- get_eic_psi_pound(Pi_star, tau_star, g, theta, psi_pound_est, S, A, Y, n)

  psi_pound_se <- sqrt(var(psi_pound_eic, na.rm = TRUE)/n)
  psi_pound_ci_lower <- psi_pound_est-1.96*psi_pound_se
  psi_pound_ci_upper <- psi_pound_est+1.96*psi_pound_se

  return(list(psi_pound_est = psi_pound_est,
              psi_pound_ci_lower = psi_pound_ci_lower,
              psi_pound_ci_upper = psi_pound_ci_upper))
}

coverage <- rep(NA, 200)
est <- rep(NA, 200)
ci_lower <- rep(NA, 200)
ci_upper <- rep(NA, 200)

for (i in 1:200) {
  print(i)
  res <- test(n, bA)
  if (res$psi_pound_ci_lower <= truth & res$psi_pound_ci_upper >= truth) {
    coverage[i] <- 1
    print("covered")
  } else {
    coverage[i] <- 0
    print("not covered")
  }
  est[i] <- res$psi_pound_est
  ci_lower[i] <- res$psi_pound_ci_lower
  ci_upper[i] <- res$psi_pound_ci_upper
}

# plot the sampling distribution
p <- ggplot(data = data.frame(est), aes(x = est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 15) +
  geom_vline(aes(xintercept = truth), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) +
  labs(title = "",
       x = "Point estimate",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

bias <- mean(est - truth)
variance <- var(est)
mse <- bias^2+variance

# save(list = c("p", "coverage", "est", "ci_lower", "ci_upper"),
#      file = "out/r_loss_psi_pound.RData")
