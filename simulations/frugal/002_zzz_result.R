library(ggplot2)
library(ggpubr)
#load("out/scm_replicate.RData")

bias_seq <- 0.5#seq(0, 2, 0.1)
ate <- 0.205

# get bias as multiples of sd of estimates
const_multi <- (apply(est_rwd, 2, mean)-ate)/apply(est_rwd, 2, sd)

# get mse
rmse_atmle_glmnet <- apply(est_atmle_glmnet, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_escvtmle <- apply(est_escvtmle, 2, function(x) sqrt(mean((x-ate)^2+var(x))))
rmse_rct_only <- apply(est_rct_only, 2, function(x) sqrt(mean((x-ate)^2+var(x))))

# get bias
bias_atmle_glmnet <- apply(est_atmle_glmnet, 2, mean)-ate
bias_escvtmle <- apply(est_escvtmle, 2, mean)-ate
bias_rct_only <- apply(est_rct_only, 2, mean)-ate

# get variance
var_atmle_glmnet <- apply(est_atmle_glmnet, 2, var)
var_escvtmle <- apply(est_escvtmle, 2, var)
var_rct_only <- apply(est_rct_only, 2, var)

# plot
par(mfrow = c(1, 3))

# 1. relative RMSE
p1 <- ggplot(data.frame(bias_seq, rmse_atmle_glmnet/rmse_rct_only, rmse_escvtmle/rmse_rct_only), aes(x = const_multi)) +
  geom_line(aes(y = rmse_atmle_glmnet/rmse_rct_only, color = "ATMLE"), linewidth = 1.5) +
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

plts <- ggarrange(p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE)
ggsave(filename = "001_scm_replicate.pdf", plot = plts, device = "pdf",
       path = "fig", width = 12, height = 4, dpi = 300)
