source("run.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
dgp_name <- "vdl25"
res_df <- run(n = 1000,
              sim_data = sim_data_vdl25,
              sim_data_args = list(bias = "b"),
              family = "gaussian",
              seed = 123,
              truth = 1.5)

saveRDS(res_df,
        file = "out/dgp_" %+% dgp_name %+% "_" %+% timestamp %+% ".RDS")
