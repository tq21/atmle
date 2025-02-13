source("utils.R")
res_df <- run(B = 20,
              bias = "c",
              g_rct = 0.67,
              controls_only = FALSE,
              ate = 4.2,
              n_rct_seq = 500,
              n_rwd_seq = 1500)

write.csv(res_df, "out/run_bias_c.csv")
