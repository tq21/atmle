library(ggplot2)

generate_data <- function(N, bA){
  # U
  UY <- rnorm(N, 0, 1)

  # W
  W1 <- rnorm(N)
  W2 <- rnorm(N)
  W3 <- rnorm(N)
  W4 <- rnorm(N)

  # A
  A <- rbinom(N, 1, 0.5)

  # Y
  Y <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+UY

  # data
  data <- data.frame(W1, W2, W3, W4, A, Y)

  return(data)
}

n <- 500
bA <- 1.3

test <- function(n) {
  data <- generate_data(n, bA)
  W_node <- c(1, 2, 3, 4)
  A_node <- 5
  Y_node <- 6
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]

  # theta_tilde(W)=E(Y|W)
  # use lasso
  theta_tilde <- learn_theta_tilde(W, Y, method = "lasso")
  g <- rep(0.5, n)
  psi <- learn_psi_tilde(W, A, Y, g, theta_tilde)
  psi_est <- mean(psi$pred)
  psi_eic <- get_eic_psi_tilde(psi, g, theta_tilde, Y, A, n)
  psi_se <- sqrt(var(psi_eic, na.rm = TRUE)/n)
  psi_ci_lower <- psi_est-1.96*psi_se
  psi_ci_upper <- psi_est+1.96*psi_se

  return(list(psi_est = psi_est,
              psi_ci_lower = psi_ci_lower,
              psi_ci_upper = psi_ci_upper))
}

coverage <- rep(NA, 200)
est <- rep(NA, 200)
ci_lower <- rep(NA, 200)
ci_upper <- rep(NA, 200)

for (i in 1:200) {
  print(i)
  res <- test(500)
  if (res$psi_ci_lower <= bA & res$psi_ci_upper >= bA) {
    coverage[i] <- 1
    print("covered")
  } else {
    coverage[i] <- 0
    print("not covered")
  }
  est[i] <- res$psi_est
  ci_lower[i] <- res$psi_ci_lower
  ci_upper[i] <- res$psi_ci_upper
}

# plot the sampling distribution
p <- ggplot(data = data.frame(est), aes(x = est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = 1.3), color = "red", linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)) +
  labs(title = "",
       x = "Point estimate",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

save(list = c("p", "coverage", "est", "ci_lower", "ci_upper"),
     file = "out/r_loss_rct.RData")
