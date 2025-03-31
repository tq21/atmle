source("sim_data.R")
res_df <- read.csv("out/res_0306_153829.csv")
truth <- get_truth()

res_df %>%
  summarize(psi_os_mse = mean((psi_os-truth)^2),
            psi_reg_cv_mse = mean((psi_reg_cv-truth)^2),
            psi_reg_mse = mean((psi_reg-truth)^2),
            psi_relax_mse = mean((psi_relax-truth)^2),
            cover_proj_reg = mean(truth >= lower_proj_reg & truth <= upper_proj_reg),
            cover_proj_relax = mean(truth >= lower_proj_relax & truth <= upper_proj_relax),
            cover_delta_reg = mean(truth >= lower_delta_reg & truth <= upper_delta_reg),
            cover_delta_relax = mean(truth >= lower_delta_relax & truth <= upper_delta_relax),
            .by = "n")

mean((res_df$psi_reg-truth)^2)
mean((res_df$psi_relax-truth)^2)
