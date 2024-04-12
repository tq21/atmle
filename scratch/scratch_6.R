sim_data <- function(ate,
                     n,
                     rct,
                     g_rct,
                     bias,
                     controls_only) {
  # error
  UY <- rnorm(n, 0, 1)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)

  # study indicator S and treatment A
  if (rct) {
    S <- rep(1, n)
    A <- rbinom(n, 1, g_rct)
  } else {
    S <- rep(0, n)
    if (controls_only) {
      A <- rep(0, n)
    } else {
      A <- rbinom(n, 1, plogis(0.5*W1))
    }
  }

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


#set.seed(123)

g_rct <- 0.5

# simulate data
data_rct <- sim_data(ate = 1.5,
                     n = 500,
                     rct = TRUE,
                     g_rct = g_rct,
                     bias = "b",
                     controls_only = FALSE)
data_rwd <- sim_data(ate = 1.5,
                     n = 500,
                     rct = FALSE,
                     g_rct = g_rct,
                     bias = "b",
                     controls_only = TRUE)
data <- rbind(data_rct, data_rwd)

S_node <- 1; W_node <- 2:4; A_node <- 5; Y_node <- 6

# apply gaussian noise to W
n_rct_A1 <- sum(data_rct$A)
W1_pert <- data_rct[data_rct$A == 1,]$W1 + rnorm(n_rct_A1, 0, 0.1)
W2_pert <- data_rct[data_rct$A == 1,]$W2 + rnorm(n_rct_A1, 0, 0.1)
W3_pert <- data_rct[data_rct$A == 1,]$W3 + rnorm(n_rct_A1, 0, 0.1)
data_perturb <- data.frame(S = 0,
                           W1 = W1_pert,
                           W2 = W2_pert,
                           W3 = W3_pert,
                           A = 1,
                           Y = data_rct[data_rct$A == 1,]$Y)
data_aug <- rbind(data, data_perturb)

theta_method <- "glm"
Pi_method <- "glm"
g_method <- "glm"
theta_tilde_method <- "glm"
Q_method <- "glm"
bias_working_model <- "glmnet"
pooled_working_model <- "glmnet"
family <- "gaussian"

# default
atmle_res <- atmle(data = data,
                   S_node = S_node,
                   W_node = W_node,
                   A_node = A_node,
                   Y_node = Y_node,
                   atmle_pooled = TRUE,
                   controls_only = TRUE,
                   theta_method = theta_method,
                   Pi_method = Pi_method,
                   g_method = g_method,
                   theta_tilde_method = theta_tilde_method,
                   Q_method = Q_method,
                   bias_working_model = bias_working_model,
                   pooled_working_model = pooled_working_model,
                   g_rct = g_rct,
                   family = family,
                   verbose = FALSE,
                   var_method = "ic")

# perturbed
atmle_res_perturb <- atmle(data = data_aug,
                           S_node = S_node,
                           W_node = W_node,
                           A_node = A_node,
                           Y_node = Y_node,
                           atmle_pooled = TRUE,
                           controls_only = FALSE,
                           theta_method = theta_method,
                           Pi_method = Pi_method,
                           g_method = g_method,
                           theta_tilde_method = theta_tilde_method,
                           Q_method = Q_method,
                           bias_working_model = bias_working_model,
                           pooled_working_model = pooled_working_model,
                           g_rct = g_rct,
                           family = family,
                           verbose = FALSE,
                           var_method = "ic")

# compare
atmle_res$upper-atmle_res$lower
atmle_res_perturb$upper-atmle_res_perturb$lower
abs(atmle_res$est-1.5)
abs(atmle_res_perturb$est-1.5)









