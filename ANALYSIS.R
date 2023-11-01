options(sl3.verbose = TRUE) # global option to print progress for sl3 fitting

atmle_fit <- atmle(data = cohort_comb,
                   S_node = S_node, # CHANGE THIS ACCORDINGLY
                   W_node = W_node, # CHANGE THIS ACCORDINGLY
                   A_node = A_node, # CHANGE THIS ACCORDINGLY
                   Y_node = Y_node, # CHANGE THIS ACCORDINGLY
                   controls_only = TRUE,
                   atmle_pooled = TRUE,
                   theta_method = "sl3",
                   Pi_method = "sl3",
                   g_method = "sl3",
                   theta_tilde_method = "sl3",
                   Q_method = "sl3",
                   bias_working_model = "glmnet",
                   pooled_working_model = "glmnet",
                   g_rct = 0.5,
                   verbose = FALSE)
print(atmle_fit)

Q.SL.library <- c("SL.speedglm", "SL.mean")
g.SL.library <- c("SL.speedglm", "SL.gam", "SL.earth")
folds <- 5

es_fit <- EScvtmle::ES.cvtmle(txinrwd = FALSE,
                              data = cohort_comb,
                              study = "STUDY",
                              covariates = covariates,
                              treatment_var = "ARM",
                              treatment = 1,
                              outcome = "MACE_ACD_evnt",
                              pRCT = 0.5,
                              V = folds,
                              Q.SL.library = Q.SL.library,
                              g.SL.library = g.SL.library,
                              Q.discreteSL = TRUE,
                              g.discreteSL = TRUE,
                              family = "binomial",
                              fluctuation = "logistic",
                              comparisons = list(c(1), c(1,0)),
                              target.gwt = TRUE,
                              cvControl = list(V = folds))
