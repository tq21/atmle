library(EScvtmle)
library(devtools)
load_all()
options(sl3.verbose = TRUE)
set.seed(2024)

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

# 1. ESCVTMLE
res_escvtmle <- ES.cvtmle(txinrwd=TRUE,
                          data=dat, study="study",
                          covariates=c("aged", "sex", "momedu", "hfiacat"),
                          treatment_var="intervention", treatment=1,
                          outcome="laz", NCO="Nlt18scale",
                          Delta=NULL, Delta_NCO=NULL,
                          pRCT=0.5, V=5, Q.SL.library=c("SL.glm"),
                          g.SL.library=c("SL.glm", "SL.earth", "SL.gam", "SL.xgboost"),
                          Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)

# 2. atmle
res_atmle_HAL <- atmle(data = dat,
                       S_node = 2,
                       W_node = c(4, 5, 6, 7),
                       A_node = 1,
                       Y_node = 3,
                       controls_only = FALSE,
                       family = "gaussian",
                       atmle_pooled = TRUE,
                       theta_method = "sl3",
                       Pi_method = "sl3",
                       g_method = "sl3",
                       theta_tilde_method = "sl3",
                       Q_method = "glm",
                       bias_working_model = "HAL",
                       pooled_working_model = "glmnet",
                       g_rct = 0.5,
                       verbose = FALSE)

res_df_unbiased <- data.frame(estimator = c("ES-CVTMLE", "A-TMLE"),
                              est = c(res_escvtmle$ATE$b2v, res_atmle_HAL$est),
                              lower = c(as.numeric(res_escvtmle$CI$b2v[1]), res_atmle_HAL$lower),
                              upper = c(as.numeric(res_escvtmle$CI$b2v[2]), res_atmle_HAL$upper))

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
                          g.SL.library=c("SL.glm", "SL.earth", "SL.gam", "SL.xgboost"),
                          Q.discreteSL=TRUE, g.discreteSL=TRUE,
                          family="gaussian", family_nco="gaussian", fluctuation = "logistic",
                          comparisons = list(c(1),c(1,0)), adjustnco = FALSE, target.gwt = TRUE)

# 2. atmle
res_atmle_HAL <- atmle(data = dat,
                       S_node = 2,
                       W_node = c(4, 5, 6, 7),
                       A_node = 1,
                       Y_node = 3,
                       controls_only = FALSE,
                       family = "gaussian",
                       atmle_pooled = TRUE,
                       theta_method = "sl3",
                       Pi_method = "sl3",
                       g_method = "sl3",
                       theta_tilde_method = "sl3",
                       Q_method = "glm",
                       bias_working_model = "HAL",
                       pooled_working_model = "glmnet",
                       g_rct = 0.5,
                       verbose = FALSE)

res_df_biased <- data.frame(estimator = c("ES-CVTMLE", "A-TMLE"),
                            est = c(res_escvtmle$ATE$b2v, res_atmle_HAL$est),
                            lower = c(as.numeric(res_escvtmle$CI$b2v[1]), res_atmle_HAL$lower),
                            upper = c(as.numeric(res_escvtmle$CI$b2v[2]), res_atmle_HAL$upper))

save(list = c("res_df_unbiased", "res_df_biased"), file = "out/wash_results.RData")
