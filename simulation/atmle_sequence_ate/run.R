.libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
            .libPaths()))
library(devtools)
library(tmle)
library(torch)
library(hal9001)
library(furrr)
library(doMC)
library(origami)
load_all()
source("sim_data.R")
set.seed(4905090)
plan(multisession, workers = availableCores()-1)
registerDoMC(cores = availableCores()-1)
ate <- get_truth()
`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

# simulation parameters
B <- 200
n_seq <- seq(500, 2000, 500)

res_df <- map_dfr(n_seq, function(n) {
  map_dfr(seq(B), function(b) {
    print("n = " %+% n %+% ", B = " %+% b)
    data <- sim_data(n)
    W <- colnames(data)[grep("W", colnames(data))]

    # SEQ-A-TMLE
    res_atmle <- atmle_ate_torch(data = data,
                                 W = W,
                                 A = "A",
                                 Y = "Y",
                                 eic_method = "svd_pseudo_inv",
                                 lr = 1e-2,
                                 family = "gaussian",
                                 browse = FALSE,
                                 parallel = TRUE)
    tmp_df <- map_dfr(res_atmle, function(.x) {
      return(data.frame(psi = .x$psi,
                        lower = .x$lower,
                        upper = .x$upper,
                        PnEIC = .x$PnEIC,
                        sn = .x$sn))
    })
    min_var_idx <- which.min(tmp_df$upper - tmp_df$lower)

    return(data.frame(b = b,
                      n = n,
                      selector = c("relax", "min_var"),
                      psi = c(tmp_df$psi[1], tmp_df$psi[min_var_idx]),
                      lower = c(tmp_df$lower[1], tmp_df$lower[min_var_idx]),
                      upper = c(tmp_df$upper[1], tmp_df$upper[min_var_idx]),
                      PnEIC = c(tmp_df$PnEIC[1], tmp_df$PnEIC[min_var_idx]),
                      sn = c(tmp_df$sn[1], tmp_df$sn[min_var_idx])))
  })
})

write.csv(res_df, "out/atmle_seq_ate" %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
