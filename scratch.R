library(devtools)
library(EScvtmle)
library(sl3)
load_all()
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

controls_only <- FALSE

data_rct <- sim_data(ate = 1.5,
                     n = 400,
                     rct = TRUE,
                     g_rct = 0.67,
                     bias = "a",
                     controls_only = controls_only)
data_rwd <- sim_data(ate = 1.5,
                     n = 2000,
                     rct = FALSE,
                     g_rct = 0.67,
                     bias = "a",
                     controls_only = controls_only)
data <- rbind(data_rct, data_rwd)

res <- atmle(data = data,
             S = "S",
             W = c("W1", "W2", "W3"),
             A = "A",
             Y = "Y",
             controls_only = controls_only,
             family = "gaussian",
             theta_method = "glm",
             g_method = "glm",
             theta_tilde_method = "glm",
             bias_working_model = "HAL",
             pooled_working_model = "HAL",
             max_degree = 1,
             verbose = FALSE)

res_escvtmle <- ES.cvtmle(txinrwd = !controls_only,
                          data = data,
                          study = "S",
                          covariates = c("W1", "W2", "W3"),
                          treatment_var = "A",
                          treatment = 1,
                          outcome = "Y",
                          pRCT = 0.67,
                          family = "gaussian",
                          Q.SL.library = c("SL.glm"),
                          g.SL.library = c("SL.glm"),
                          Q.discreteSL = TRUE,
                          g.discreteSL = TRUE,
                          V = 5)

res$upper-res$lower
as.numeric(res_escvtmle$CI$b2v[2]-res_escvtmle$CI$b2v[1])
