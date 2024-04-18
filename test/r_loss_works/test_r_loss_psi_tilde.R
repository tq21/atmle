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
    #A <- rbinom(N, 1, plogis(-0.1+0.2*W1+0.5*W2-0.1*W3))
    A <- rbinom(N, 1, 0.5)
  }

  # S
  if (is.null(S)) {
    S <- rbinom(N, 1, 0.2)
  }

  # Y
  Y <- 0.3+bA*A+0.5*W1+0.3*W3-0.5*W4+UY

  # data
  data <- data.frame(S, W1, W2, W3, W4, A, Y)

  return(data)
}

n <- 1000
bA <- 1.3
# truth <- mean(generate_data(100000, bA, NULL, rep(1, 100000))$Y) -
#   mean(generate_data(100000, bA, NULL, rep(0, 100000))$Y)
#truth <- bA + 10
truth <- 1.3

test <- function(n) {
  data <- generate_data(n, bA)
  S_node <- 1
  W_node <- c(2, 3, 4, 5)
  A_node <- 6
  Y_node <- 7
  S <- data[, S_node]
  W <- data[, W_node]
  A <- data[, A_node]
  Y <- data[, Y_node]

  # learn nuisance
  theta_tilde <- learn_theta_tilde(W, Y, method = "glm")
  g <- learn_g(S, W, A, 0.5, method = "glm")

  # learn psi_tilde using R-loss
  psi_tilde <- learn_psi_tilde(W, A, Y, g, theta_tilde, method = "lasso")
  psi_tilde_est <- mean(psi_tilde$pred)
  psi_tilde_eic <- get_eic_psi_tilde(psi_tilde, g, theta_tilde, Y, A, n)

  psi_tilde_se <- sqrt(var(psi_tilde_eic, na.rm = TRUE)/n)
  psi_tilde_ci_lower <- psi_tilde_est-1.96*psi_tilde_se
  psi_tilde_ci_upper <- psi_tilde_est+1.96*psi_tilde_se

  return(list(est = psi_tilde_est,
              ci_lower = psi_tilde_ci_lower,
              ci_upper = psi_tilde_ci_upper))
}

coverage <- rep(NA, 200)
est <- rep(NA, 200)
ci_lower <- rep(NA, 200)
ci_upper <- rep(NA, 200)

for (i in 1:200) {
  print(i)
  res <- test(n)
  if (res$ci_lower <= truth & res$ci_upper >= truth) {
    coverage[i] <- 1
    print("covered")
  } else {
    coverage[i] <- 0
    print("not covered")
  }
  est[i] <- res$est
  ci_lower[i] <- res$ci_lower
  ci_upper[i] <- res$ci_upper
}

# plot the sampling distribution
p <- ggplot(data = data.frame(est), aes(x = est)) +
  geom_histogram(fill = "#0072B2", color = "white", bins = 20) +
  geom_vline(aes(xintercept = truth), color = "red", linetype = "dashed", linewidth = 1) +
  #scale_x_continuous(breaks = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6)) +
  labs(title = "",
       x = "Point estimate",
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12))

save(list = c("p", "coverage", "est", "ci_lower", "ci_upper"),
     file = "out/r_loss_psi_pound.RData")
