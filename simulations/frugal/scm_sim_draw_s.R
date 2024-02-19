library(atmle)
library(EScvtmle)
library(devtools)
load_all()

`%+%` <- function(a, b) paste0(a, b)

# DGP: SCM VERSION
sim_data <- function(n, gamma, rct_prop) {
  W1 <- rnorm(n, 0, 1)
  W2 <- rbinom(n, 1, 0.5)
  W3 <- rnorm(n, 0.1*W2, 1)
  U1 <- rnorm(n, 0, 1)
  U2 <- rnorm(n, 0, 0.35)

  S <- rbinom(n, 1, rct_prop)

  A <- numeric(length = n)
  A[S == 1] <- rbinom(sum(S), 1, 0.5)
  A[S == 0] <- rbinom(n - sum(S), 1, plogis(W1+W2+gamma*U1+gamma*U2))

  # Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)
  Y <- rnorm(n, 0.2*A+0.1*W3+U1+U2, 1)

  data <- data.frame(S = S,
                     W1 = W1,
                     W2 = W2,
                     W3 = W3,
                     A = A,
                     Y = Y)

  return(data)
}

# SCM SIMULATIONS
set.seed(124)
B <- 100
n <- 10000
bias_seq <- c(0.5)#seq(0, 5, 0.1)
g_rct <- 0.5
rct_prop <- 0.5
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6
est_atmle_glmnet <- matrix(nrow = B, ncol = length(bias_seq))
est_escvtmle <- matrix(nrow = B, ncol = length(bias_seq))
est_rct_only <- matrix(nrow = B, ncol = length(bias_seq))
est_rwd <- matrix(nrow = B, ncol = length(bias_seq))

for (b in 1:B) {
  print("run: " %+% b)

  # simulate data
  data_list <- lapply(bias_seq, function(bias) {
    return(sim_data(n, bias, rct_prop))
  })

  for (i in 1:length(bias_seq)) {
    data <- data_list[[i]]

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
                         Qform = "Y ~ A + A*W3",
                         g.SL.library = c("SL.glm"),
                         Q.SL.library = c("SL.glm"),
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
                              covariates = c("W1", "W2", "W3"),
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
              "est_rwd"), file = "scm_sim_draw_s.RData")
