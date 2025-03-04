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
data <- sim_data(500)
W <- colnames(data)[grep("W", colnames(data))]

# SEQ-A-TMLE
set.seed(123)
res_svd_pseudo_inv <- atmle_ate_torch(data = data,
                                      W = W,
                                      A = "A",
                                      Y = "Y",
                                      eic_method = "svd_pseudo_inv",
                                      lr = 1e-2,
                                      family = "gaussian",
                                      browse = FALSE,
                                      parallel = TRUE)
set.seed(123)
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
  data.frame(eic_method = c("svd-based pseudoinverse", "add 1e-3 to diagonal"),
             psi = c(.x$psi, .y$psi),
             lower = c(.x$lower, .y$lower),
             upper = c(.x$upper, .y$upper),
             PnEIC = c(.x$PnEIC, .y$PnEIC),
             sn = c(.x$sn, .y$sn))
})
idx <- which.min(res_seq_atmle_df$upper - res_seq_atmle_df$lower)
res_seq_atmle_df$index <- seq_along(res_seq_atmle_df$psi)

fit_1 <- ggplot(res_seq_atmle_df, aes(x = index, y = psi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~eic_method) +
  geom_hline(yintercept = truth, linetype = "dashed", color = "red") +
  labs(x = "lambda index", y = "Psi Estimate",
       title = "") +
  theme_bw()

fig_2 <- ggplot(res_seq_atmle_df, aes(x = index)) +
  geom_line(aes(y = abs(PnEIC), color = "|PnEIC|")) +
  geom_line(aes(y = sn, color = "threshold")) +
  facet_wrap(~eic_method) +
  labs(x = "lambda index", y = "", color = "") +
  scale_color_manual(values = c("|PnEIC|" = "black", "threshold" = "blue")) +
  theme_minimal() +
  theme_bw()

fig <- ggarrange(fit_1, fig_2, common.legend = TRUE, legend = "bottom", ncol = 1)

ggsave(filename = "diag_vs_svd_ate.pdf", plot = fig, device = "pdf",
       path = "figs", width = 5, height = 5, dpi = 300)
