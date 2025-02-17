.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(furrr)
library(doMC)
library(origami)
load_all()
source("sim_data.R")
set.seed(923845)
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)
ate <- get_truth()
`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

# simulation parameters
B <- 200
n_seq <- seq(500, 2000, 500)

res_df <- map_dfr(n_seq, function(n) {
  map_dfr(seq(B), function(b) {
    print("n = " %+% n %+% ", B = " %+% b)
    data <- sim_data(n)
    W <- colnames(data)[grep("W", colnames(data))]

    # SEQ-A-TMLE
    res_atmle <- atmle_ate_torch(data = data,
                                 W = W,
                                 A = "A",
                                 Y = "Y",
                                 eic_method = "svd_pseudo_inv",
                                 lr = 1e-2,
                                 family = "gaussian",
                                 browse = FALSE,
                                 parallel = TRUE)
    tmp_df <- map_dfr(res_atmle, function(.x) {
      return(data.frame(psi = .x$psi,
                        lower = .x$lower,
                        upper = .x$upper,
                        PnEIC = .x$PnEIC,
                        sn = .x$sn))
    })
    min_var_idx <- which.min(tmp_df$upper - tmp_df$lower)

    # HAL-TMLE
    g_basis_list <- enumerate_basis(x = data[, W, drop = FALSE],
                                    max_degree = 2,
                                    smoothness_orders = 1)
    g_design_mat <- make_design_matrix(X = as.matrix(data[, W, drop = FALSE]),
                                       blist = g_basis_list)
    g_fit <- cv.glmnet(x = g_design_mat,
                       y = data$A,
                       family = "binomial",
                       alpha = 1,
                       nlambda = 50,
                       nfolds = 10,
                       parallel = TRUE)
    g1W <- as.numeric(predict(g_fit, newx = g_design_mat, s = "lambda.min", type = "response"))
    Q_basis_list <- enumerate_basis(x = data[, c(W, "A"), drop = FALSE],
                                    max_degree = 2,
                                    smoothness_orders = 1)
    Q1_counter <- cbind(data[, W, drop = FALSE], A = 1)
    Q0_counter <- cbind(data[, W, drop = FALSE], A = 0)
    Q_design_mat <- make_design_matrix(X = as.matrix(data[, c(W, "A"), drop = FALSE]),
                                       blist = Q_basis_list)
    Q1_design_mat <- make_design_matrix(X = as.matrix(Q1_counter),
                                       blist = Q_basis_list)
    Q0_design_mat <- make_design_matrix(X = as.matrix(Q0_counter),
                                       blist = Q_basis_list)
    Q_fit <- cv.glmnet(x = Q_design_mat,
                       y = data$Y,
                       family = "gaussian",
                       alpha = 1,
                       nlambda = 50,
                       nfolds = 10,
                       parallel = TRUE)
    Q0_init <- as.numeric(predict(Q_fit, newx = Q0_design_mat, s = "lambda.min"))
    Q1_init <- as.numeric(predict(Q_fit, newx = Q1_design_mat, s = "lambda.min"))
    tmle_res <- tmle(W = as.matrix(data[, W, drop = FALSE]),
                     A = as.numeric(data$A),
                     Y = as.numeric(data$Y),
                     family = "gaussian",
                     Q = cbind(Q0_init, Q1_init),
                     g1W = g1W)

    return(data.frame(b = b,
                      n = n,
                      selector = c("relax", "min_var", "TMLE"),
                      psi = c(tmp_df$psi[1], tmp_df$psi[min_var_idx], tmle_res$estimates$ATE$psi),
                      lower = c(tmp_df$lower[1], tmp_df$lower[min_var_idx], tmle_res$estimates$ATE$CI[1]),
                      upper = c(tmp_df$upper[1], tmp_df$upper[min_var_idx], tmle_res$estimates$ATE$CI[2]),
                      PnEIC = c(tmp_df$PnEIC[1], tmp_df$PnEIC[min_var_idx], mean(tmle_res$estimates$IC$IC.ATE)),
                      sn = c(tmp_df$sn[1], tmp_df$sn[min_var_idx], NA)))
  })
})

write.csv(res_df, "out/atmle_seq_ate" %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
