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
truth_tilde <- mean(generate_data(100000, bA, NULL, rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, NULL, rep(0, 100000))$Y)
S1_true_ate_A <- mean(generate_data(100000, bA, rep(1, 100000), rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, rep(1, 100000), rep(0, 100000))$Y)
S0_true_ate_A <- mean(generate_data(100000, bA, rep(0, 100000), rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, rep(0, 100000), rep(0, 100000))$Y)
truth_pound <- 0.5*S0_true_ate_A-0.5*S1_true_ate_A
truth <- truth_tilde - truth_pound
truth <- bA + 10

test <- function(n, bA) {
  data <- generate_data(n, bA)
  res <- atmle(data,
               S_node = 1,
               W_node = c(2, 3, 4, 5),
               A_node = 6,
               Y_node = 7,
               verbose = FALSE)
  return(res)
}

coverage <- rep(NA, 200)
est <- rep(NA, 200)
ci_lower <- rep(NA, 200)
ci_upper <- rep(NA, 200)

tilde_coverage <- rep(NA, 200)
tilde_est <- rep(NA, 200)
tilde_ci_lower <- rep(NA, 200)
tilde_ci_upper <- rep(NA, 200)

pound_coverage <- rep(NA, 200)
pound_est <- rep(NA, 200)
pound_ci_lower <- rep(NA, 200)
pound_ci_upper <- rep(NA, 200)

for (i in 1:200) {
  print(i)
  res <- test(n, bA)
  if (res$lower <= truth & res$upper >= truth) {
    coverage[i] <- 1
    print("psi covered")
  } else {
    coverage[i] <- 0
    print("psi not covered")
  }

  if (res$psi_pound_lower <= truth_pound & res$psi_pound_upper >= truth_pound) {
    pound_coverage[i] <- 1
    print("pound covered")
  } else {
    pound_coverage[i] <- 0
    print("pound not covered")
  }

  if (res$psi_tilde_lower <= truth_tilde & res$psi_tilde_upper >= truth_tilde) {
    tilde_coverage[i] <- 1
    print("tilde covered")
  } else {
    tilde_coverage[i] <- 0
    print("tilde not covered")
  }

  est[i] <- res$est
  ci_lower[i] <- res$lower
  ci_upper[i] <- res$upper

  pound_est[i] <- res$psi_pound_est
  pound_ci_lower[i] <- res$psi_pound_lower
  pound_ci_upper[i] <- res$psi_pound_upper

  tilde_est[i] <- res$psi_tilde_est
  tilde_ci_lower[i] <- res$psi_tilde_lower
  tilde_ci_upper[i] <- res$psi_tilde_upper
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
#      file = "out/atmle_both_constant_bias.RData")
