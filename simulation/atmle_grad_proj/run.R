`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
library(devtools)
library(dplyr)
library(furrr)
library(doMC)
load_all()
source("sim_data.R")
set.seed(12345)
registerDoMC(cores = availableCores()-1)

B <- 100
n_seq <- seq(500, 2000, 500)
res_df <- map_dfr(n_seq, function(.n) {
  map_dfr(seq(B), function(.b) {
    print("n: " %+% .n %+% ", run: " %+% .b)
    data <- sim_data(.n)
    res <- atmle_single(data = data,
                        W = c("W1", "W2", "W3"),
                        A = "A",
                        Y = "Y",
                        family = "gaussian",
                        theta_method = "glm",
                        g_method = "glm",
                        tau_A_method = "HAL",
                        enumerate_basis_args = list(max_degree = 2,
                                                    smoothness_orders = 1),
                        browse = FALSE)

    return(data.frame(n = .n,
                      b = .b,
                      psi_reg = res$psi_reg,
                      psi_relax = res$psi_relax,
                      lower_proj_reg = res$lower_proj_reg,
                      lower_proj_relax = res$lower_proj_relax,
                      upper_proj_reg = res$upper_proj_reg,
                      upper_proj_relax = res$upper_proj_relax,
                      lower_delta_reg = res$lower_delta_reg,
                      lower_delta_relax = res$lower_delta_relax,
                      upper_delta_reg = res$upper_delta_reg,
                      upper_delta_relax = res$upper_delta_relax))
  })
})

write.csv(res_df, file = "out/res_" %+% timestamp %+% ".csv", row.names = FALSE)
