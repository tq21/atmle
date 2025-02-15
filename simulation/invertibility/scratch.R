library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(EScvtmle)
library(furrr)
library(doMC)
library(ggplot2)
load_all()
source("sim_data.R")
#set.seed(4905090)
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)

truth <- get_truth()

set.seed(982475)
data <- sim_data(500)
W <- colnames(data)[grep("W", colnames(data))]

# SEQ-A-TMLE
res_svd_pseudo_inv <- atmle_ate_torch(data = data,
                                      W = W,
                                      A = "A",
                                      Y = "Y",
                                      eic_method = "svd_pseudo_inv",
                                      lr = 1e-2,
                                      family = "gaussian",
                                      browse = FALSE,
                                      parallel = TRUE)
res_diag <- atmle_ate_torch(data = data,
                            W = W,
                            A = "A",
                            Y = "Y",
                            eic_method = "diag",
                            lr = 1e-2,
                            family = "gaussian",
                            browse = FALSE,
                            parallel = TRUE)


res_seq_atmle_df <- map2_dfr(res_svd_pseudo_inv, res_diag, function(.x, .y) {
  data.frame(eic_method = c("svd_pseudo_inv", "diag"),
             psi = c(.x$psi, .y$psi),
             lower = c(.x$lower, .y$lower),
             upper = c(.x$upper, .y$upper),
             PnEIC = c(.x$PnEIC, .y$PnEIC),
             sn = c(.x$sn, .y$sn))
})
idx <- which.min(res_seq_atmle_df$upper - res_seq_atmle_df$lower)
res_seq_atmle_df$index <- seq_along(res_seq_atmle_df$psi)

ggplot(res_seq_atmle_df, aes(x = index, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~eic_method) +
  geom_hline(yintercept = truth, linetype = "dashed", color = "red") +
  labs(x = "Index", y = "Psi Estimate",
       title = "Sequence of Psi Estimates with 95% Confidence Intervals") +
  theme_minimal()

abs(res_seq_atmle_df$PnEIC) <= res_seq_atmle_df$sn






