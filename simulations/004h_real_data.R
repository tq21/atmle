library(ggpubr)

options(sl3.verbose = TRUE)
source("utils.R")

set.seed(2022)

# RCT: study 1
# unbiased external: study 2
# biased external: study 3

# unbiased external
data(wash)
dat <- wash[which(wash$study %in% c(1, 2)),]
dat$study[which(dat$study == 2)] <- 0
dat$sex <- as.numeric(dat$sex) - 1
dat$momedu <- as.numeric(dat$momedu) - 1
dat$hfiacat <- as.numeric(dat$hfiacat) - 1

# 1. ESCVTMLE (benchmark)
res_escvtmle <- ES.cvtmle(txinrwd=TRUE,
                          data=dat, study="study",
                          covariates=c("aged", "sex", "momedu", "hfiacat"),
                          treatment_var="intervention", treatment=1,
                          outcome="laz", NCO="Nlt18scale",
                          Delta=NULL, Delta_NCO=NULL,
                          pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
                          g.SL.library=c("SL.glm", "SL.earth"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)

# 2. atmle (with main term lasso working model)
res_atmle_glmnet <- atmle(data = dat,
                          S_node = 2,
                          W_node = c(4, 5, 6, 7),
                          A_node = 1,
                          Y_node = 3,
                          controls_only = FALSE,
                          family = "gaussian",
                          atmle_pooled = TRUE,
                          theta_method = "glm",
                          Pi_method = "glm",
                          g_method = "glm",
                          theta_tilde_method = "glm",
                          Q_method = "glm",
                          bias_working_model = "glmnet",
                          pooled_working_model = "glmnet",
                          g_rct = 0.5,
                          verbose = TRUE)

res_df_unbiased <- data.frame(estimator = c("ES-CVTMLE", "A-TMLE"),
                              est = c(res_escvtmle$ATE$b2v, res_atmle_glmnet$est),
                              lower = c(as.numeric(res_escvtmle$CI$b2v[1]), res_atmle_glmnet$lower),
                              upper = c(as.numeric(res_escvtmle$CI$b2v[2]), res_atmle_glmnet$upper))

# biased external
data(wash)
dat <- wash[which(wash$study %in% c(1, 3)),]
dat$study[which(dat$study == 3)] <- 0
dat$sex <- as.numeric(dat$sex) - 1
dat$momedu <- as.numeric(dat$momedu) - 1
dat$hfiacat <- as.numeric(dat$hfiacat) - 1

# 1. ESCVTMLE (benchmark)
res_escvtmle <- ES.cvtmle(txinrwd=TRUE,
                          data=dat, study="study",
                          covariates=c("aged", "sex", "momedu", "hfiacat"),
                          treatment_var="intervention", treatment=1,
                          outcome="laz", NCO="Nlt18scale",
                          Delta=NULL, Delta_NCO=NULL,
                          pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
                          g.SL.library=c("SL.glm", "SL.earth"), Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)

# 2. atmle (with main term lasso working model)
res_atmle_glmnet <- atmle(data = dat,
                          S_node = 2,
                          W_node = c(4, 5, 6, 7),
                          A_node = 1,
                          Y_node = 3,
                          controls_only = FALSE,
                          family = "gaussian",
                          atmle_pooled = TRUE,
                          theta_method = "glm",
                          Pi_method = "glm",
                          g_method = "glm",
                          theta_tilde_method = "glm",
                          Q_method = "glm",
                          bias_working_model = "glmnet",
                          pooled_working_model = "glmnet",
                          g_rct = 0.5,
                          verbose = TRUE)

res_df_biased <- data.frame(estimator = c("ES-CVTMLE", "A-TMLE"),
                            est = c(res_escvtmle$ATE$b2v, res_atmle_glmnet$est),
                            lower = c(as.numeric(res_escvtmle$CI$b2v[1]), res_atmle_glmnet$lower),
                            upper = c(as.numeric(res_escvtmle$CI$b2v[2]), res_atmle_glmnet$upper))

save(list = c("res_df_unbiased", "res_df_biased"), file = "out/wash.RData")

res_df_unbiased$bias <- "unbiased external"
res_df_biased$bias <- "biased external"
res_df <- rbind(res_df_unbiased, res_df_biased)

# plot (unbiased external)
p_unbiased <- ggplot(res_df_unbiased, aes(x = estimator, y = est)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "dodgerblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5, color = "dodgerblue") +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.1),
                     labels = c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "Unbiased external", x = "", y = "Estimate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16))

# plot (biased external)
p_biased <- ggplot(res_df_biased, aes(x = estimator, y = est)) +
  geom_point(position = position_dodge(width = 0.5), size = 3.5, color = "darkorange") +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.5), width = 0.3, size = 1.5, color = "darkorange") +
  scale_y_continuous(limits = c(-0.3, 0.3), breaks = seq(-0.3, 0.3, 0.1),
                     labels = c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.5) +
  labs(title = "Biased external", x = "", y = "Estimate", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16))

wash_plt <- ggarrange(p_unbiased, p_biased, nrow = 1, ncol = 2)
ggsave(filename = "wash.pdf", plot = wash_plt, device = "pdf",
       path = "plot", width = 10, height = 8, dpi = 300)
