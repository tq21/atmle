library(ggplot2)

load_all()

generate_data <- function(N, bA, S=NULL, A=NULL){
  # U
  UY <- rnorm(N, 0, 1)

  # S
  if (is.null(S)) {
    S <- rbinom(N, 1, 0.2)
  }

  # W
  W1 <- vector(length = N)
  W1[S == 0] <- rnorm(N - sum(S), mean = 1, sd = 1)
  W1[S == 1] <- rnorm(sum(S), mean = 1, sd = 1)

  W2 <- vector(length = N)
  W2[S == 0] <- rnorm(N - sum(S), mean = 1, sd = 1)
  W2[S == 1] <- rnorm(sum(S), mean = 1, sd = 1)

  W3 <- vector(length = N)
  W3[S == 0] <- rnorm(N - sum(S), mean = 1, sd = 1)
  W3[S == 1] <- rnorm(sum(S), mean = 1, sd = 1)

  W4 <- vector(length = N)
  W4[S == 0] <- rnorm(N - sum(S), mean = 1, sd = 1)
  W4[S == 1] <- rnorm(sum(S), mean = 1, sd = 1)

  # A
  if (is.null(A)) {
    A <- vector(length = N)
    A[S == 0] <- rbinom(N - sum(S), 1, plogis(-0.1+0.2*W1+0.5*W2-0.1*W3))
    A[S == 1] <- rbinom(sum(S), 1, 0.5)
  }

  # Y
  Y <- vector(length = N)
  Y_S0 <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+0.7+UY # bias = 0.7
  Y_S1 <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+UY
  Y[S == 0] <- Y_S0[S == 0]
  Y[S == 1] <- Y_S1[S == 1]

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

# simulation parameters --------------------------------------------------------
B <- 200
n <- 1000
bA <- 1.3
truth <- bA
truth_tilde <- mean(generate_data(100000, bA, NULL, rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, NULL, rep(0, 100000))$Y)
truth_pound <- mean(generate_data(100000, bA, rep(1, 100000), rep(1, 100000))$Y) -
  mean(generate_data(100000, bA, rep(1, 100000), rep(0, 100000))$Y)
# truth <- mean(generate_data(100000, bA, NULL, rep(1, 100000))$Y) -
#   mean(generate_data(100000, bA, NULL, rep(0, 100000))$Y)

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

# psi
psi_coverage <- rep(NA, B)
psi_est <- rep(NA, B)
psi_ci_lower <- rep(NA, B)
psi_ci_upper <- rep(NA, B)

# psi tilde
tilde_coverage <- rep(NA, B)
tilde_est <- rep(NA, B)
tilde_ci_lower <- rep(NA, B)
tilde_ci_upper <- rep(NA, B)

# psi pound
pound_coverage <- rep(NA, B)
pound_est <- rep(NA, B)
pound_ci_lower <- rep(NA, B)
pound_ci_upper <- rep(NA, B)

for (i in 1:B) {
  print(i)
  res <- test(n, bA)
  if (res$lower <= truth & res$upper >= truth) {
    psi_coverage[i] <- 1
    print("psi covered")
  } else {
    psi_coverage[i] <- 0
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

  psi_est[i] <- res$est
  psi_ci_lower[i] <- res$lower
  psi_ci_upper[i] <- res$upper

  pound_est[i] <- res$psi_pound_est
  pound_ci_lower[i] <- res$psi_pound_lower
  pound_ci_upper[i] <- res$psi_pound_upper

  tilde_est[i] <- res$psi_tilde_est
  tilde_ci_lower[i] <- res$psi_tilde_lower
  tilde_ci_upper[i] <- res$psi_tilde_upper
}

# plot the sampling distributions
p_psi <- ggplot(data = data.frame(psi_est), aes(x = psi_est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = truth), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)) +
  labs(title = "",
       x = "Point estimate (Psi)",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

p_tilde <- ggplot(data = data.frame(tilde_est), aes(x = tilde_est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = truth_tilde), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)) +
  labs(title = "",
       x = "Point estimate (Psi_tilde)",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

p_pound <- ggplot(data = data.frame(pound_est), aes(x = pound_est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = truth_pound), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)) +
  labs(title = "",
       x = "Point estimate (Psi_pound)",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

save(list = c("psi_coverage", "psi_est", "psi_ci_lower", "psi_ci_upper",
              "tilde_coverage", "tilde_est", "tilde_ci_lower", "tilde_ci_upper",
              "pound_coverage", "pound_est", "pound_ci_lower", "pound_ci_upper",
              "p_psi", "p_tilde", "p_pound", "generate_data"),
     file = "out/atmle_lasso_SWAY_no_bias.RData")
