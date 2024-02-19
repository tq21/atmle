library(ggplot2)
library(ggpubr)
load("out/001_replicate_HAL.RData")

bias_seq <- seq(0, 5, 0.2)
ate <- 0.205

# get bias as multiples of sd of estimates
const_multi <- abs((apply(est_pooled, 2, mean)-ate))/apply(est_pooled, 2, sd)

# get mse
rmse_atmle <- apply(est_atmle, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_escvtmle <- apply(est_escvtmle, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_rct_only <- apply(est_rct_only, 2, function(x) sqrt(mean((x-ate)^2+var(x))))

# get bias
bias_atmle <- apply(est_atmle, 2, mean)-ate
bias_escvtmle <- apply(est_escvtmle, 2, mean)-ate
bias_rct_only <- apply(est_rct_only, 2, mean)-ate

# get variance
var_atmle <- apply(est_atmle, 2, var)
var_escvtmle <- apply(est_escvtmle, 2, var)
var_rct_only <- apply(est_rct_only, 2, var)

# get bias of the bias
psi_pound <- apply(psi_pound, 2, mean)
atmle_psi_pound <- apply(atmle_psi_pound, 2, mean)
bias_of_bias <- (psi_pound-atmle_psi_pound)^2

# plot
par(mfrow = c(1, 3))

# 1. relative RMSE
p1 <- ggplot(data.frame(bias_seq, rmse_atmle/rmse_rct_only, rmse_escvtmle/rmse_rct_only), aes(x = const_multi)) +
  geom_line(aes(y = rmse_atmle/rmse_rct_only, color = "ATMLE"), linewidth = 1.5) +
  geom_line(aes(y = rmse_escvtmle/rmse_rct_only, color = "ES-CVTMLE"), linewidth = 1.5) +
  geom_line(aes(y = 1, color = "RCT only"), linewidth = 1.5) +
  scale_color_manual(name = "Estimator", values = c("red", "orange", "black")) +
  labs(x = "Bias as a multiple of sd(ATE)", y = "Relative RMSE") +
  theme_bw()

p2 <- ggplot(data.frame(bias_seq, bias_atmle_glmnet^2, bias_escvtmle^2, bias_rct_only^2), aes(x = const_multi)) +
  geom_line(aes(y = bias_atmle_glmnet^2, color = "ATMLE"), linewidth = 1.5) +
  geom_line(aes(y = bias_escvtmle^2, color = "ES-CVTMLE"), linewidth = 1.5) +
  geom_line(aes(y = bias_rct_only^2, color = "RCT only"), linewidth = 1.5) +
  scale_color_manual(name = "Estimator", values = c("red", "orange", "black")) +
  labs(x = "Bias as a multiple of sd(ATE)", y = "Bias^2") +
  theme_bw()

p3 <- ggplot(data.frame(bias_seq, var_atmle_glmnet, var_escvtmle, var_rct_only), aes(x = const_multi)) +
  geom_line(aes(y = var_atmle_glmnet, color = "ATMLE"), linewidth = 1.5) +
  geom_line(aes(y = var_escvtmle, color = "ES-CVTMLE"), linewidth = 1.5) +
  geom_line(aes(y = var_rct_only, color = "RCT only"), linewidth = 1.5) +
  scale_color_manual(name = "Estimator", values = c("red", "orange", "black")) +
  labs(x = "Bias as a multiple of sd(ATE)", y = "Variance") +
  theme_bw()

p4 <- ggplot(data.frame(bias_seq, bias_of_bias), aes(x = const_multi)) +
  geom_line(aes(y = bias_of_bias, color = "Bias of bias"), linewidth = 1.5) +
  scale_color_manual(name = "Estimator", values = c("red", "orange", "black")) +
  labs(x = "Bias as a multiple of sd(ATE)", y = "Bias of bias") +
  theme_bw()

plts <- ggarrange(p1, p2, p3, p4, nrow = 1, ncol = 4, common.legend = TRUE)
# ggsave(filename = "scm.pdf", plot = plts, device = "pdf",
#        path = "fig", width = 20, height = 4, dpi = 300)
