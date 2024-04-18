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

  Y <- rnorm(n, (0.2+0.1*W3)*A+U1+U2, 1)

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
n <- 1500
bias_seq <- seq(0, 2, 0.1)
g_rct <- 0.5
rct_prop <- 0.2
controls_only <- FALSE
S_node <- 1
W_node <- 2:4
A_node <- 5
Y_node <- 6
atmle_psi_pound <- matrix(nrow = B, ncol = length(bias_seq))
psi_pound <- matrix(nrow = B, ncol = length(bias_seq))

for (b in 1:B) {
  print("run: " %+% b)

  # simulate data
  data_list <- lapply(bias_seq, function(bias) {
    return(sim_data(n, bias, rct_prop))
  })

  for (i in 1:length(bias_seq)) {
    data <- data_list[[i]]

    # nonparametric
    # nonparam_res <- nonparametric(data = data,
    #                               S_node = S_node,
    #                               W_node = W_node,
    #                               A_node = A_node,
    #                               Y_node = Y_node,
    #                               controls_only = controls_only,
    #                               family = "gaussian",
    #                               g_rct = g_rct,
    #                               verbose = FALSE)
    nonparam_res <- rct_only(data = data,
                             S_node = S_node,
                             W_node = W_node,
                             A_node = A_node,
                             Y_node = Y_node,
                             g_rct = g_rct,
                             nuisance_method = "glm",
                             family = "gaussian",
                             verbose = FALSE)

    # tmle on pooled-data
    tmle_pooled_res <- tmle(Y = data[, Y_node],
                            A = data[, A_node],
                            W = data[, W_node],
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

    # collect results
    atmle_psi_pound[b,i] <- atmle_glmnet_res$psi_pound_est
    psi_pound[b,i] <- tmle_pooled_res$estimates$ATE$psi-nonparam_res$est
  }
}

save(list = c("atmle_psi_pound",
              "psi_pound",), file = "out/psi_pound_draw_s.RData")
