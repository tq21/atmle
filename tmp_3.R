library(devtools)
load_all()
sim_data <- function(ate,
                     n,
                     rct,
                     g,
                     bias) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)

  # study indicator S
  if (rct) S <- rep(1, n) else S <- rep(0, n)

  # treatment A
  A <- rbinom(n, 1, g)

  # bias term for RWD data
  if (bias == "a") {
    b <- 0.2+0.1*W1*(1-A)
  } else if (bias == "b") {
    b <- 0.5+3.1*W1*(1-A)+0.8*W3
  }

  # outcome
  Y <- 2.5+0.9*W1+1.1*W2+2.7*W3+ate*A+UY+(1-S)*b

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

data_rct <- sim_data(ate = 1.5,
                     n = 500,
                     rct = TRUE,
                     g = 0.5,
                     bias = "a")
data_rwd <- sim_data(ate = 1.5,
                     n = 4000,
                     rct = FALSE,
                     g = 0.5,
                     bias = "a")
data <- rbind(data_rct, data_rwd)
res_robust <- atmle(data = data,
                    robust_covariate = TRUE,
                    S = "S",
                    W = c("W1", "W2", "W3"),
                    A = "A",
                    Y = "Y",
                    controls_only = FALSE,
                    family = "gaussian",
                    theta_method = "glm",
                    Pi_method = "glm",
                    g_method = "glm",
                    theta_tilde_method = "glm",
                    bias_working_model = "glmnet",
                    pooled_working_model = "glmnet",
                    verbose = FALSE)
res <- atmle(data = data,
             robust_covariate = FALSE,
             S = "S",
             W = c("W1", "W2", "W3"),
             A = "A",
             Y = "Y",
             controls_only = FALSE,
             family = "gaussian",
             theta_method = "glm",
             Pi_method = "glm",
             g_method = "glm",
             theta_tilde_method = "glm",
             bias_working_model = "glmnet",
             pooled_working_model = "glmnet",
             verbose = FALSE)

res_robust$upper-res_robust$lower
res$upper-res$lower
