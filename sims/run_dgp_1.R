source("run.R")
source("dgps/sim_data.R")
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
res_df <- run(n = 1000,
              sim_data = sim_data,
              sim_data_args = list(bias = "b"),
              family = "gaussian",
              seed = 123,
              truth = 1.5,
              truth_2 = 1.655082)

saveRDS(res_df,
        file = "out/_" %+% timestamp %+% ".RDS")
