library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(EScvtmle)
library(furrr)
library(doMC)
library(ggplot2)
load_all()
source("sim_data_fusion.R")
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)
set.seed(4905090)

# simulate data
ate <- 1.5
g_rct <- 0.67
controls_only <- FALSE
data_rct <- sim_data(ate = ate,
                     n = 400,
                     rct = TRUE,
                     g_rct = g_rct,
                     bias = "b",
                     controls_only = controls_only)
data_rwd <- sim_data(ate = ate,
                     n = 400,
                     rct = FALSE,
                     g_rct = g_rct,
                     bias = "b",
                     controls_only = controls_only)
data <- rbind(data_rct, data_rwd)

# A-TMLE, diag method
set.seed(123)
res_diag <- atmle_torch(data = data,
                        S = "S",
                        W = c("W1", "W2", "W3"),
                        A = "A",
                        Y = "Y",
                        controls_only = controls_only,
                        family = "gaussian",
                        eic_method = "diag",
                        parallel = TRUE,
                        browse = FALSE)

# A-TMLE, svd method
set.seed(123)
res_svd_pseudo_inv <- atmle_torch(data = data,
                                  S = "S",
                                  W = c("W1", "W2", "W3"),
                                  A = "A",
                                  Y = "Y",
                                  controls_only = controls_only,
                                  family = "gaussian",
                                  eic_method = "svd_pseudo_inv",
                                  parallel = TRUE,
                                  browse = FALSE)

# ES-CVTMLE
res_escvtmle <- ES.cvtmle(txinrwd = !controls_only,
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
res_df <- map2_dfr(res_svd_pseudo_inv, res_diag, function(.x, .y) {
  data.frame(est_name = c("svd_pseudo_inv", "diag"),
             psi = c(.x$psi, .y$psi),
             lower = c(.x$lower, .y$lower),
             upper = c(.x$upper, .y$upper),
             PnEIC = c(.x$PnEIC, .y$PnEIC),
             sn = c(.x$sn, .y$sn))
})
res_df <- rbind(res_df, data.frame(est_name = "ES-CVTMLE",
                                   psi = res_escvtmle$ATE$b2v,
                                   lower = res_escvtmle$CI$b2v[1],
                                   upper = res_escvtmle$CI$b2v[2],
                                   PnEIC = NA,
                                   sn = NA))
res_df$index <- seq_along(res_df$psi)

# plot
ggplot(res_df, aes(x = index, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~est_name) +
  geom_hline(yintercept = ate, linetype = "dashed", color = "red") +
  labs(x = "Index", y = "Estimated Psi",
       title = "") +
  theme_minimal()
