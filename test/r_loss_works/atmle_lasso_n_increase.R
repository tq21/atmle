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
  Y <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+UY # no bias

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

# simulation parameters --------------------------------------------------------
B <- 200
n <- 500
bA <- 1.5
truth <- bA

test <- function(n, bA) {
  data <- generate_data(n, bA)
  res <- atmle(data,
               S_node = 1,
               W_node = c(2, 3, 4, 5),
               A_node = 6,
               Y_node = 7,
               nuisance_method="lasso",
               working_model="HAL",
               verbose = FALSE)
  return(res)
}

# psi
psi_coverage <- rep(NA, B)
psi_est <- rep(NA, B)
psi_ci_lower <- rep(NA, B)
psi_ci_upper <- rep(NA, B)

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

  psi_est[i] <- res$est
  psi_ci_lower[i] <- res$lower
  psi_ci_upper[i] <- res$upper
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

save(list = c("psi_coverage", "psi_est", "psi_ci_lower", "psi_ci_upper",
              "p_psi", "generate_data"),
     file = "out/atmle_hal_sway_no_bias.RData")
