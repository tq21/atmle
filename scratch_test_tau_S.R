library(devtools)
library(EScvtmle)
library(sl3)
load_all()
#set.seed(123)
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

set.seed(123)
controls_only <- FALSE

data_rct <- sim_data(ate = 1.5,
                     n = 400,
                     rct = TRUE,
                     g_rct = 0.67,
                     bias = "b",
                     controls_only = controls_only)
data_rwd <- sim_data(ate = 1.5,
                     n = 2000,
                     rct = FALSE,
                     g_rct = 0.67,
                     bias = "b",
                     controls_only = controls_only)
data <- rbind(data_rct, data_rwd)
S <- data$S
W <- data[, grep("W", colnames(data))]
A <- data$A
Y <- data$Y
folds <- make_folds(n = length(S), V = 10)
foldid <- unlist(map(folds, function(.fold) {
  rep(.fold$v, length(.fold$validation_set))
}))
idx <- unlist(map(folds, function(.fold) {
  .fold$validation_set
}))
foldid <- foldid[idx]

g1W <- learn_g(S = S,
               W = W,
               A = A,
               method = "glm",
               controls_only = controls_only,
               folds = folds,
               g_bounds = c(0, 1),
               cross_fit_nuisance = TRUE)
Pi <- learn_Pi(g = g1W,
               A = A,
               Pi_bounds = c(0, 1))
theta_WA <- learn_theta_WA(W = W,
                           A = A,
                           Y = Y,
                           delta = rep(1, length(Y)),
                           controls_only = controls_only,
                           method = "glm",
                           folds = folds,
                           family = "gaussian",
                           theta_bounds = NULL,
                           cross_fit_nuisance = TRUE)

tau <- learn_tau_S(S = S,
                   W = W,
                   A = A,
                   Y = Y,
                   Pi = Pi,
                   theta = theta_WA,
                   g1W = g1W,
                   delta = rep(1, length(Y)),
                   controls_only = controls_only,
                   method = "HAL",
                   v_folds = 5,
                   min_working_model = NULL,
                   target_gwt = FALSE,
                   Pi_bounds = c(0, 1),
                   enumerate_basis_args = list(max_degree = 2,
                                               smoothness_orders = 1,
                                               num_knots = c(20, 5)),
                   weights = rep(1, length(Y)),
                   bias_working_model_formula = NULL,
                   target_method = "oneshot",
                   verbose = TRUE)
