source("utils.R")

meta <- list(ate = 1.5,
             g_rct = 0.5,
             bias = "small_bias",
             gamma = 1,
             nuisance_method = "glm",
             working_model = "glmnet",
             B = 1000,
             sim_fun = run_sim,
             sim_data = sim_data)

results <- run_sim(ate = meta$ate,
                   g_rct = meta$g_rct,
                   bias = meta$bias,
                   gamma = meta$gamma,
                   nuisance_method = meta$nuisance_method,
                   working_model = meta$working_model,
                   B = meta$B)

save(list = c("meta", "results"),
     file = "out/1120_gamma_1.RData")
