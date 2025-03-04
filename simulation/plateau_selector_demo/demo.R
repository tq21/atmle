library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(EScvtmle)
library(furrr)
library(doMC)
library(ggplot2)
library(ggpubr)
load_all()
source("sim_data.R")
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)

set.seed(4905090)
truth <- get_truth()
data <- sim_data(1000)
W <- colnames(data)[grep("W", colnames(data))]

# SEQ-A-TMLE
set.seed(123)
res <- atmle_ate_torch(data = data,
                       W = W,
                       A = "A",
                       Y = "Y",
                       eic_method = "svd_pseudo_inv",
                       lr = 1e-2,
                       family = "gaussian",
                       browse = FALSE,
                       parallel = TRUE)

res_df <- map_dfr(res, function(.x) {
  data.frame(psi = .x$psi,
             var = var(.x$EIC),
             lower = .x$lower,
             upper = .x$upper,
             PnEIC = .x$PnEIC,
             L1_norm = .x$L1_norm,
             sn = .x$sn)
})
res_df$index <- seq_along(res_df$psi)

# enforce monotonicity
basis_list <- enumerate_basis(x = res_df$L1_norm,
                              smoothness_orders = 0)
design_mat <- make_design_matrix(X = as.matrix(res_df$L1_norm),
                                 blist = basis_list)
fit_mono <- cv.glmnet(x = design_mat,
                      y = res_df$var,
                      family = "gaussian",
                      alpha = 1,
                      lower.limits = 0,
                      parallel = TRUE)
pred_mono <- predict(fit_mono, s = "lambda.min", newx = design_mat)

