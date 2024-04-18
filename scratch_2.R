set.seed(1234)

# nuisance: glm; working model: glmnet
atmle_res <- atmle(data = as.data.frame(cohort_comb),
                   S_node = which(names(cohort_comb) == "STUDY"), #"STUDY"
                   W_node = which(names(cohort_comb) %in% covariates),
                   A_node = which(names(cohort_comb) == "ARM"), # "ARM"
                   Y_node = which(names(cohort_comb) == "MACE_ACD_evnt"),
                   atmle_pooled = TRUE,
                   controls_only = T,
                   theta_method = "glm",
                   Pi_method = "glm",
                   g_method = "glm",
                   theta_tilde_method = "glm",
                   bias_working_model = "glmnet",
                   pooled_working_model = "glmnet",
                   g_rct = 0.5,
                   family = "binomial",
                   v_folds = folds,
                   verbose = T)

# nuisance: glm; working model: glmnet
# with censored outcomes, please encode missing Y as NA in the data
atmle_res <- atmle(data = as.data.frame(cohort_comb),
                   S_node = which(names(cohort_comb) == "STUDY"), #"STUDY"
                   W_node = which(names(cohort_comb) %in% covariates),
                   A_node = which(names(cohort_comb) == "ARM"), # "ARM"
                   Y_node = which(names(cohort_comb) == "MACE_ACD_evnt"),
                   atmle_pooled = TRUE,
                   controls_only = T,
                   theta_method = "glm",
                   Pi_method = "glm",
                   g_method = "glm",
                   g_delta_method = "glm",
                   theta_tilde_method = "glm",
                   bias_working_model = "glmnet",
                   pooled_working_model = "glmnet",
                   g_rct = 0.5,
                   family = "binomial",
                   v_folds = folds,
                   verbose = T)

# nuisance: sl3; working model: glmnet
atmle_res <- atmle(data = as.data.frame(cohort_comb),
                   S_node = which(names(cohort_comb) == "STUDY"), #"STUDY"
                   W_node = which(names(cohort_comb) %in% covariates),
                   A_node = which(names(cohort_comb) == "ARM"), # "ARM"
                   Y_node = which(names(cohort_comb) == "MACE_ACD_evnt"),
                   atmle_pooled = TRUE,
                   controls_only = T,
                   theta_method = "sl3",
                   Pi_method = "sl3",
                   g_method = "sl3",
                   theta_tilde_method = "sl3",
                   bias_working_model = "glmnet",
                   pooled_working_model = "glmnet",
                   g_rct = 0.5,
                   family = "binomial",
                   v_folds = folds,
                   verbose = T)

# nuisance: sl3; working model: glmnet
# with censored outcomes, please encode missing Y as NA in the data
atmle_res <- atmle(data = as.data.frame(cohort_comb),
                   S_node = which(names(cohort_comb) == "STUDY"), #"STUDY"
                   W_node = which(names(cohort_comb) %in% covariates),
                   A_node = which(names(cohort_comb) == "ARM"), # "ARM"
                   Y_node = which(names(cohort_comb) == "MACE_ACD_evnt"),
                   atmle_pooled = TRUE,
                   controls_only = T,
                   theta_method = "sl3",
                   Pi_method = "sl3",
                   g_method = "sl3",
                   g_delta_method = "sl3",
                   theta_tilde_method = "glm",
                   bias_working_model = "glmnet",
                   pooled_working_model = "glmnet",
                   g_rct = 0.5,
                   family = "binomial",
                   v_folds = folds,
                   verbose = T)
