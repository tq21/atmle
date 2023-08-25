generate_data <- function(N, bA){
  # U
  UY <- rnorm(N, 0, 1)

  # W
  W1 <- rnorm(N)
  W2 <- rnorm(N)
  W3 <- rnorm(N)
  W4 <- rnorm(N)

  # A
  A <- rbinom(N, 1, plogis(-0.1+0.2*W1+0.5*W2-0.1*W3))

  # S
  S <- rbinom(N, 1, plogis(0.2-0.5*W1-0.3*W2+0.2*W3-0.1*A))

  # Y
  Y <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

n <- 1000
bA <- 1.3
bS <- 0

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

for (i in 1:200) {
  print(i)
  res <- test(n, bA)
  if (res$lower <= bA & res$upper >= bA) {
    coverage[i] <- 1
    print("covered")
  } else {
    coverage[i] <- 0
    print("not covered")
  }
  est[i] <- res$est
  ci_lower[i] <- res$lower
  ci_upper[i] <- res$upper
}

# plot the sampling distribution
p <- ggplot(data = data.frame(est), aes(x = est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = bA), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2)) +
  labs(title = "",
       x = "Point estimate",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

save(list = c("p", "coverage", "est", "ci_lower", "ci_upper"),
     file = "out/r_loss_combined_no_bias_psi_pound.RData")
