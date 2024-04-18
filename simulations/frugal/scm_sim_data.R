library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

sim_data <- function(ate, n, rct_prop, g_rct, bias, controls_only) {
  # error
  UY <- rnorm(n, 0, 1)
  Ubias <- rnorm(n, 0, 0.1)

  # baseline covariates
  W1 <- rnorm(n, 0, 1)
  W2 <- rnorm(n, 0, 1)
  W3 <- rnorm(n, 0, 1)
  W4 <- rnorm(n, 0, 1)

  # study indicator, S=1 for RCT, S=0 for RWD
  S <- rbinom(n, 1, rct_prop)

  # treatments (external data has both treated and controls)
  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, g_rct)
  if (controls_only) {
    A[S == 0] <- rep(0, n - sum(S))
  } else {
    A[S == 0] <- rbinom(n - sum(S), 1, plogis(0.8*W1+0.9*W2))
  }

  # outcome
  Y <- -2.4-1.9*W1-2.5*W2+1.2*W3+2.7*W4+ate*A+UY+(1-S)*bias*A+(1-S)*Ubias

  # data frames combining RCT and RWD
  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     W4 = W4,
                     A = A,
                     Y = Y)

  return(data)
}


# SCM SIMULATIONS
set.seed(124)
B <- 100
ate <- 1.5
bias_seq <- seq(0, 5, 0.1)
rct_prop <- 0.2
g_rct <- 0.5
n <- 2000
controls_only <- FALSE
S_node <- 1
W_node <- 2:5
A_node <- 6
Y_node <- 7
est_atmle_glmnet <- matrix(nrow = B, ncol = length(bias_seq))
est_escvtmle <- matrix(nrow = B, ncol = length(bias_seq))
est_rct_only <- matrix(nrow = B, ncol = length(bias_seq))
est_rwd <- matrix(nrow = B, ncol = length(bias_seq))

for (b in 1:B) {
  print("run: " %+% b)

  # simulate data
  dat.list <- lapply(bias_seq, function(bias) {
    return(sim_data(ate, n, rct_prop, g_rct, bias, controls_only))
  })

  for (i in 1:length(bias_seq)) {
    data <- as.data.frame(dat.list[[i]])

    # rct only
    rct_only_res <- rct_only(data = data,
                             S_node = S_node,
                             W_node = W_node,
                             A_node = A_node,
                             Y_node = Y_node,
                             g_rct = g_rct,
                             nuisance_method = "glm",
                             family = "gaussian",
                             verbose = FALSE)

    # tmle on RWD
    tmle_rwd_res <- tmle(Y = data[data$S == 0, Y_node],
                         A = data[data$S == 0, A_node],
                         W = data[data$S == 0, W_node],
                         g.SL.library = c("SL.glm", "SL.earth", "SL.gam"),
                         Q.SL.library = c("SL.glm", "SL.earth", "SL.gam"),
                         family = "gaussian")

    # atmle with glmnet
    atmle_glmnet_res <- atmle(data = data,
                              S_node = S_node,
                              W_node = W_node,
                              A_node = A_node,
                              Y_node = Y_node,
                              atmle_pooled = TRUE,
                              controls_only = FALSE,
                              theta_method = "glm",
                              Pi_method = "glm",
                              g_method = "glm",
                              theta_tilde_method = "glm",
                              Q_method = "glm",
                              bias_working_model = "glmnet",
                              pooled_working_model = "glmnet",
                              g_rct = g_rct,
                              family = "gaussian",
                              verbose = FALSE)

    # ES-CVTMLE
    escvtmle_res <- ES.cvtmle(txinrwd = TRUE,
                              data = data,
                              study = "S",
                              covariates = c("W1", "W2", "W3", "W4"),
                              treatment_var = "A",
                              treatment = 1,
                              outcome = "Y",
                              pRCT = g_rct,
                              family = "gaussian",
                              Q.SL.library = c("SL.glm"),
                              g.SL.library = c("SL.glm"),
                              Q.discreteSL = TRUE,
                              g.discreteSL = TRUE,
                              V = 5)

    # collect results
    est_atmle_glmnet[b,i] <- atmle_glmnet_res$est
    est_rct_only[b,i] <- rct_only_res$est
    est_escvtmle[b,i] <- escvtmle_res$ATE$b2v
    est_rwd[b,i] <- tmle_rwd_res$estimates$ATE$psi
  }
}

save(list = c("est_atmle_glmnet",
              "est_escvtmle",
              "est_rct_only",
              "est_rwd"), file = "scm_sim_results.RData")
